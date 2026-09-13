//! Debug and terrain rendering (design doc section 23).
//!
//! The point of these layers is diagnosis. When a world looks wrong, the
//! question is always *which stage* is wrong - plate motion, boundary
//! detection, classification, deformation, or erosion - and each of those has
//! a layer here that shows it directly rather than through its effect on the
//! final terrain.

use crate::boundary::BoundaryKind;
use crate::math::{Vec2, percentile, smoothstep};
use crate::simulation::Simulation;
use crate::tectonics::falloff;
use image::{Rgb, RgbImage};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layer {
    /// Finished map: hypsometric tint, hillshade, rivers and lakes.
    Terrain,
    /// Elevation alone, no water and no shading.
    Elevation,
    /// Plate ownership, with sites and velocity arrows.
    Plates,
    /// Boundary classification: red convergent, blue divergent, green transform.
    Boundaries,
    /// Where deformation is actually being applied, and how hard.
    Strength,
    /// Crust age, young to old.
    CrustAge,
    /// Rainfall.
    Rainfall,
    /// Flow accumulation, on a log scale.
    Flow,
}

impl Layer {
    pub const ALL: [Layer; 8] = [
        Layer::Terrain,
        Layer::Elevation,
        Layer::Plates,
        Layer::Boundaries,
        Layer::Strength,
        Layer::CrustAge,
        Layer::Rainfall,
        Layer::Flow,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Layer::Terrain => "terrain",
            Layer::Elevation => "elevation",
            Layer::Plates => "plates",
            Layer::Boundaries => "boundaries",
            Layer::Strength => "strength",
            Layer::CrustAge => "crust-age",
            Layer::Rainfall => "rainfall",
            Layer::Flow => "flow",
        }
    }

    pub fn parse(s: &str) -> Option<Layer> {
        Layer::ALL
            .into_iter()
            .find(|l| l.name() == s || l.name().replace('-', "_") == s)
    }
}

pub fn render(sim: &Simulation, layer: Layer) -> RgbImage {
    match layer {
        Layer::Terrain => terrain(sim, true),
        Layer::Elevation => terrain(sim, false),
        Layer::Plates => plates(sim),
        Layer::Boundaries => boundaries(sim),
        Layer::Strength => strength(sim),
        Layer::CrustAge => crust_age(sim),
        Layer::Rainfall => rainfall(sim),
        Layer::Flow => flow(sim),
    }
}

// --- layers ---------------------------------------------------------------

/// Hypsometric terrain. With `water`, hillshade, rivers and lakes go on top.
fn terrain(sim: &Simulation, water: bool) -> RgbImage {
    let w = &sim.world;
    let sea = sim.sea_level();
    let mut img = RgbImage::new(w.width as u32, w.height as u32);

    // Rivers are drawn where flow is large relative to the map's own scale,
    // so the threshold does not need retuning per world size.
    let river_threshold = percentile(&w.flow, 0.995).max(1.0);

    for y in 0..w.height {
        for x in 0..w.width {
            let idx = w.idx(x, y);
            let e = w.elevation[idx];
            let mut c = HYPSOMETRIC.sample(e);

            if water {
                if e > sea {
                    c = shade(c, hillshade(sim, x, y));
                }
                if w.water[idx] > 0.02 && e > sea {
                    c = mix(c, [58, 102, 160], 0.85); // lake
                }
                let f = w.flow[idx];
                if e > sea && f > river_threshold {
                    let t = smoothstep(river_threshold, river_threshold * 12.0, f);
                    c = mix(c, [46, 88, 150], 0.45 + 0.45 * t);
                }
            }
            img.put_pixel(x as u32, y as u32, Rgb(c));
        }
    }
    img
}

/// Plate ownership, plus each plate's site and velocity arrow.
fn plates(sim: &Simulation) -> RgbImage {
    let w = &sim.world;
    let mut img = RgbImage::new(w.width as u32, w.height as u32);

    for idx in 0..w.len() {
        let (x, y) = w.coords(idx);
        let id = w.plate_id[idx];
        let mut c = plate_color(id);
        // Continental plates are drawn lighter, so the crust layout reads at
        // a glance.
        if sim.plates[id as usize].crust_type.is_continental() {
            c = mix(c, [255, 255, 255], 0.35);
        }
        img.put_pixel(x as u32, y as u32, Rgb(c));
    }

    // Boundaries in black, so plate outlines stay visible over the fill.
    for b in &sim.boundaries {
        let (x, y) = w.coords(b.cell as usize);
        img.put_pixel(x as u32, y as u32, Rgb([20, 20, 20]));
    }

    // Arrows are scaled so a typical plate speed reads clearly at map size;
    // the length is proportional to speed, so relative motion stays legible.
    for p in &sim.plates {
        let tip = p.position + p.velocity * 220.0;
        draw_arrow(&mut img, p.position, tip, [15, 15, 15], 1.6);
        draw_disc(&mut img, p.position, 3.0, [250, 250, 250]);
        draw_disc(&mut img, p.position, 1.6, [20, 20, 20]);
    }
    img
}

/// Boundary classification. This is the layer that answers "is the simulation
/// deciding the right *kind* of boundary here?".
fn boundaries(sim: &Simulation) -> RgbImage {
    let w = &sim.world;
    let mut img = RgbImage::new(w.width as u32, w.height as u32);
    for idx in 0..w.len() {
        let (x, y) = w.coords(idx);
        let e = w.elevation[idx];
        // Dim terrain underneath for context.
        let base = mix(HYPSOMETRIC.sample(e), [0, 0, 0], 0.6);
        img.put_pixel(x as u32, y as u32, Rgb(base));
    }

    let peak = sim
        .boundaries
        .iter()
        .map(|b| b.strength)
        .fold(1e-6f32, f32::max);

    for b in &sim.boundaries {
        let (x, y) = w.coords(b.cell as usize);
        let hue = match b.kind {
            BoundaryKind::Convergent => [235, 70, 60],
            BoundaryKind::Divergent => [70, 130, 240],
            BoundaryKind::Transform => [90, 210, 110],
        };
        // Brightness carries strength, so a fast boundary stands out from a
        // barely-moving one of the same kind.
        let t = 0.35 + 0.65 * (b.strength / peak).clamp(0.0, 1.0);
        img.put_pixel(x as u32, y as u32, Rgb(mix([0, 0, 0], hue, t)));
    }
    img
}

/// How much deformation each cell is receiving: strength through the same
/// falloff the tectonics stage uses. Bright means "this is being worked on".
fn strength(sim: &Simulation) -> RgbImage {
    let w = &sim.world;
    let mut img = RgbImage::new(w.width as u32, w.height as u32);
    let radius = sim.field.max_radius;

    let mut influence = vec![0.0f32; w.len()];
    for idx in 0..w.len() {
        if let Some((src, d)) = sim.field.at(idx) {
            let b = &sim.boundaries[src as usize];
            influence[idx] = b.strength * falloff(d, radius);
        }
    }
    let peak = influence.iter().cloned().fold(1e-6, f32::max);

    for idx in 0..w.len() {
        let (x, y) = w.coords(idx);
        let t = (influence[idx] / peak).clamp(0.0, 1.0);
        img.put_pixel(x as u32, y as u32, Rgb(INFERNO.sample(t)));
    }
    img
}

fn crust_age(sim: &Simulation) -> RgbImage {
    let w = &sim.world;
    let oldest = percentile(&w.crust_age, 0.99).max(1.0);
    scalar_layer(sim, &w.crust_age, 0.0, oldest, &AGE)
}

fn rainfall(sim: &Simulation) -> RgbImage {
    let w = &sim.world;
    let wettest = percentile(&w.rainfall, 0.99).max(1e-3);
    scalar_layer(sim, &w.rainfall, 0.0, wettest, &RAIN)
}

/// Flow accumulation on a log scale - without the log, one trunk river
/// saturates and every tributary disappears.
fn flow(sim: &Simulation) -> RgbImage {
    let w = &sim.world;
    let logged: Vec<f32> = w.flow.iter().map(|f| (1.0 + f.max(0.0)).ln()).collect();
    let peak = percentile(&logged, 0.999).max(1e-3);
    let mut img = scalar_layer(sim, &logged, 0.0, peak, &FLOW);
    // Mask out the sea: accumulated rain in the ocean is not a river.
    for idx in 0..w.len() {
        if w.elevation[idx] <= sim.sea_level() {
            let (x, y) = w.coords(idx);
            img.put_pixel(x as u32, y as u32, Rgb([12, 16, 26]));
        }
    }
    img
}

fn scalar_layer(sim: &Simulation, data: &[f32], lo: f32, hi: f32, ramp: &ColorRamp) -> RgbImage {
    let w = &sim.world;
    let mut img = RgbImage::new(w.width as u32, w.height as u32);
    let span = (hi - lo).max(1e-6);
    for idx in 0..w.len() {
        let (x, y) = w.coords(idx);
        let t = ((data[idx] - lo) / span).clamp(0.0, 1.0);
        img.put_pixel(x as u32, y as u32, Rgb(ramp.sample(t)));
    }
    img
}

// --- shading --------------------------------------------------------------

/// Lambert shading from a north-west sun, on the elevation gradient.
fn hillshade(sim: &Simulation, x: usize, y: usize) -> f32 {
    let w = &sim.world;
    let at = |dx: i64, dy: i64| -> f32 {
        w.neighbor(x, y, dx, dy)
            .map(|i| w.elevation[i])
            .unwrap_or_else(|| w.elevation[w.idx(x, y)])
    };
    const EXAGGERATION: f32 = 2.2;
    let dzdx = (at(1, 0) - at(-1, 0)) * 0.5 * EXAGGERATION;
    let dzdy = (at(0, 1) - at(0, -1)) * 0.5 * EXAGGERATION;

    // Normal (-dzdx, -dzdy, 1) against a light from the north-west.
    let light = [-0.55f32, -0.55, 0.63];
    let len = (dzdx * dzdx + dzdy * dzdy + 1.0).sqrt();
    let dot = (-dzdx * light[0] - dzdy * light[1] + light[2]) / len;
    (0.55 + 0.75 * dot).clamp(0.35, 1.35)
}

fn shade(c: [u8; 3], factor: f32) -> [u8; 3] {
    [
        (c[0] as f32 * factor).clamp(0.0, 255.0) as u8,
        (c[1] as f32 * factor).clamp(0.0, 255.0) as u8,
        (c[2] as f32 * factor).clamp(0.0, 255.0) as u8,
    ]
}

fn mix(a: [u8; 3], b: [u8; 3], t: f32) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    [
        (a[0] as f32 + (b[0] as f32 - a[0] as f32) * t) as u8,
        (a[1] as f32 + (b[1] as f32 - a[1] as f32) * t) as u8,
        (a[2] as f32 + (b[2] as f32 - a[2] as f32) * t) as u8,
    ]
}

/// Distinct, evenly spread plate colours from the golden-ratio hue sequence.
fn plate_color(id: u16) -> [u8; 3] {
    let hue = (id as f32 * 0.618_034).fract();
    let sat = 0.45 + 0.2 * ((id as f32 * 0.37).fract());
    let val = 0.55 + 0.3 * ((id as f32 * 0.71).fract());
    hsv(hue, sat, val)
}

fn hsv(h: f32, s: f32, v: f32) -> [u8; 3] {
    let i = (h * 6.0).floor();
    let f = h * 6.0 - i;
    let p = v * (1.0 - s);
    let q = v * (1.0 - f * s);
    let t = v * (1.0 - (1.0 - f) * s);
    let (r, g, b) = match (i as i32).rem_euclid(6) {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    };
    [
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8,
    ]
}

// --- drawing --------------------------------------------------------------

fn draw_arrow(img: &mut RgbImage, from: Vec2, to: Vec2, color: [u8; 3], width: f32) {
    draw_line(img, from, to, color, width);
    let dir = (to - from).normalize_or_zero();
    if dir == Vec2::ZERO {
        return;
    }
    let head = (to - from).length().clamp(4.0, 9.0);
    let back = to - dir * head;
    let side = Vec2::new(-dir.y, dir.x) * (head * 0.45);
    draw_line(img, to, back + side, color, width);
    draw_line(img, to, back - side, color, width);
}

/// Straight line, wrapping in x and clipped at the walls.
fn draw_line(img: &mut RgbImage, from: Vec2, to: Vec2, color: [u8; 3], width: f32) {
    let steps = ((to - from).length().ceil() as i32).max(1);
    for i in 0..=steps {
        let t = i as f32 / steps as f32;
        let p = from + (to - from) * t;
        if width <= 1.0 {
            put_wrapped(img, p.x, p.y, color);
        } else {
            draw_disc(img, p, width * 0.5, color);
        }
    }
}

fn draw_disc(img: &mut RgbImage, centre: Vec2, radius: f32, color: [u8; 3]) {
    let r = radius.ceil() as i32;
    for dy in -r..=r {
        for dx in -r..=r {
            if ((dx * dx + dy * dy) as f32).sqrt() <= radius {
                put_wrapped(img, centre.x + dx as f32, centre.y + dy as f32, color);
            }
        }
    }
}

fn put_wrapped(img: &mut RgbImage, x: f32, y: f32, color: [u8; 3]) {
    let w = img.width() as i64;
    let h = img.height() as i64;
    let yi = y.round() as i64;
    if yi < 0 || yi >= h {
        return;
    }
    let xi = (x.round() as i64).rem_euclid(w);
    img.put_pixel(xi as u32, yi as u32, Rgb(color));
}

// --- colour ramps ---------------------------------------------------------

pub struct ColorRamp {
    stops: &'static [(f32, [u8; 3])],
}

impl ColorRamp {
    /// Colour at `v`, interpolating between the surrounding stops.
    pub fn sample(&self, v: f32) -> [u8; 3] {
        let stops = self.stops;
        if v <= stops[0].0 {
            return stops[0].1;
        }
        for pair in stops.windows(2) {
            let (v0, c0) = pair[0];
            let (v1, c1) = pair[1];
            if v <= v1 {
                let t = if (v1 - v0).abs() < 1e-9 {
                    0.0
                } else {
                    (v - v0) / (v1 - v0)
                };
                return mix(c0, c1, t);
            }
        }
        stops[stops.len() - 1].1
    }
}

/// Elevation tint. Stops are in the simulation's elevation units, with sea
/// level at 0.
pub const HYPSOMETRIC: ColorRamp = ColorRamp {
    stops: &[
        (-9.0, [8, 24, 58]),
        (-5.0, [14, 52, 104]),
        (-1.5, [32, 96, 158]),
        (-0.25, [86, 148, 194]),
        (0.0, [206, 194, 148]),
        (0.35, [104, 148, 84]),
        (1.6, [142, 158, 92]),
        (3.0, [152, 130, 88]),
        (4.8, [124, 106, 90]),
        (6.4, [170, 162, 156]),
        (7.6, [246, 246, 248]),
    ],
};

const INFERNO: ColorRamp = ColorRamp {
    stops: &[
        (0.0, [8, 8, 20]),
        (0.25, [72, 20, 92]),
        (0.5, [176, 48, 78]),
        (0.75, [238, 122, 32]),
        (1.0, [252, 252, 190]),
    ],
};

const AGE: ColorRamp = ColorRamp {
    stops: &[
        (0.0, [250, 248, 168]),
        (0.35, [228, 128, 70]),
        (0.7, [150, 52, 96]),
        (1.0, [34, 24, 72]),
    ],
};

const RAIN: ColorRamp = ColorRamp {
    stops: &[
        (0.0, [246, 238, 208]),
        (0.3, [162, 206, 156]),
        (0.6, [58, 150, 168]),
        (1.0, [24, 52, 128]),
    ],
};

const FLOW: ColorRamp = ColorRamp {
    stops: &[
        (0.0, [12, 16, 26]),
        (0.4, [30, 70, 110]),
        (0.75, [86, 162, 206]),
        (1.0, [232, 248, 255]),
    ],
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plate::PlateParams;
    use crate::simulation::SimulationParams;

    fn sim() -> Simulation {
        let mut s = Simulation::new(SimulationParams {
            width: 96,
            height: 48,
            seed: 5,
            plates: PlateParams {
                count: 8,
                ..Default::default()
            },
            ..Default::default()
        });
        s.run(12);
        s
    }

    #[test]
    fn every_layer_renders_at_the_world_size() {
        let sim = sim();
        for layer in Layer::ALL {
            let img = render(&sim, layer);
            assert_eq!(img.width(), 96, "{} is the wrong width", layer.name());
            assert_eq!(img.height(), 48, "{} is the wrong height", layer.name());
        }
    }

    #[test]
    fn layers_are_not_blank() {
        let sim = sim();
        for layer in Layer::ALL {
            let img = render(&sim, layer);
            let first = img.get_pixel(0, 0);
            assert!(
                img.pixels().any(|p| p != first),
                "{} rendered a flat image",
                layer.name()
            );
        }
    }

    #[test]
    fn layer_names_round_trip() {
        for layer in Layer::ALL {
            assert_eq!(Layer::parse(layer.name()), Some(layer));
        }
        assert_eq!(Layer::parse("crust_age"), Some(Layer::CrustAge));
        assert_eq!(Layer::parse("nonsense"), None);
    }

    #[test]
    fn ramps_are_clamped_at_both_ends() {
        assert_eq!(HYPSOMETRIC.sample(-500.0), [8, 24, 58]);
        assert_eq!(HYPSOMETRIC.sample(500.0), [246, 246, 248]);
        assert_eq!(INFERNO.sample(-1.0), [8, 8, 20]);
        assert_eq!(INFERNO.sample(2.0), [252, 252, 190]);
    }

    #[test]
    fn drawing_wraps_across_the_seam_instead_of_panicking() {
        let mut img = RgbImage::new(32, 16);
        draw_arrow(&mut img, Vec2::new(31.0, 8.0), Vec2::new(40.0, 8.0), [255, 0, 0], 1.0);
        draw_disc(&mut img, Vec2::new(0.0, 0.0), 3.0, [0, 255, 0]);
        assert_eq!(img.get_pixel(8, 8), &Rgb([255, 0, 0]));
        assert_eq!(img.get_pixel(30, 0), &Rgb([0, 255, 0]));
    }

}
