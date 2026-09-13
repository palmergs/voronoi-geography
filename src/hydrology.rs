//! Rainfall and water routing (design doc section 17).
//!
//! Rivers are never drawn. They emerge from flow accumulation over the
//! heightfield, which is what lets tectonics and erosion interact: a new
//! mountain range immediately changes where the water goes.
//!
//! The routing uses priority-flood with an epsilon tilt (Barnes et al.). One
//! pass fills every depression up to its outlet *and* yields a strictly
//! descending drainage direction for every cell, with no separate sink-filling
//! iteration. The order cells are popped in is kept: it is a topological order
//! of the drainage network, so flow accumulation is a single reverse sweep.

use crate::math::{Vec2, smoothstep};
use crate::noise::CylinderNoise;
use crate::world::{OFFSETS8, World};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Clone, Copy, Debug)]
pub struct HydrologyParams {
    pub sea_level: f32,
    /// Rain that falls regardless of latitude.
    pub base_rainfall: f32,
    /// Extra rain where wind climbs terrain, and the rain shadow behind it.
    pub orographic_gain: f32,
    pub rain_shadow: f32,
    /// Amplitude and scale of the noise that breaks up the latitude bands.
    pub rain_noise: f32,
    pub rain_feature_size: f32,
    /// Tilt applied inside filled depressions so water still has somewhere to go.
    pub fill_epsilon: f32,
}

impl Default for HydrologyParams {
    fn default() -> Self {
        HydrologyParams {
            sea_level: 0.0,
            base_rainfall: 0.15,
            orographic_gain: 1.6,
            rain_shadow: 1.1,
            rain_noise: 0.35,
            rain_feature_size: 90.0,
            fill_epsilon: 1e-4,
        }
    }
}

/// The drainage network for the current heightfield.
pub struct FlowNetwork {
    /// Elevation with depressions filled, tilted so every cell drains.
    pub filled: Vec<f32>,
    /// Downstream neighbour, or -1 for an outlet (ocean or map edge).
    pub receiver: Vec<i32>,
    /// Cells in ascending `filled` order - a topological order of the network.
    pub order: Vec<u32>,
}

impl FlowNetwork {
    pub fn new(n: usize) -> Self {
        FlowNetwork {
            filled: vec![0.0; n],
            receiver: vec![-1; n],
            order: Vec::with_capacity(n),
        }
    }
}

pub struct Hydrology {
    pub params: HydrologyParams,
    rain_noise: CylinderNoise,
}

impl Hydrology {
    pub fn new(params: HydrologyParams, seed: u64, width: usize) -> Self {
        Hydrology {
            rain_noise: CylinderNoise::new(
                (seed as u32) ^ 0x_5A10_BEE5,
                width,
                params.rain_feature_size,
                4,
            ),
            params,
        }
    }

    /// Fill `world.rainfall`.
    ///
    /// Three ingredients: a latitude profile (wet equator, dry subtropics, wet
    /// mid-latitudes), large-scale noise so the bands are not stripes, and an
    /// orographic term that wrings rain out of air climbing terrain and starves
    /// the far side. Prevailing wind direction comes from the latitude band,
    /// which is why mountain ranges end up with a wet and a dry flank.
    pub fn rainfall(&self, world: &mut World) {
        let p = &self.params;

        for y in 0..world.height {
            let lat = world.latitude(y);
            let band = latitude_rainfall(lat);
            let wind = prevailing_wind(lat);

            for x in 0..world.width {
                let idx = world.idx(x, y);
                let n = self.rain_noise.get(x as f32, y as f32);
                let mut rain = (band + p.rain_noise * n).max(0.0) + p.base_rainfall;

                // Air arrives from upwind; compare terrain there to here.
                if world.elevation[idx] > p.sea_level {
                    let up = world
                        .neighbor(x, y, -wind.x.signum() as i64, 0)
                        .unwrap_or(idx);
                    let climb = world.elevation[idx] - world.elevation[up].max(p.sea_level);
                    if climb > 0.0 {
                        rain *= 1.0 + p.orographic_gain * smoothstep(0.0, 2.5, climb);
                    } else {
                        rain /= 1.0 + p.rain_shadow * smoothstep(0.0, 2.5, -climb);
                    }
                }

                world.rainfall[idx] = rain;
            }
        }
    }

    /// Fill depressions, pick drainage directions, and accumulate flow.
    ///
    /// Writes `world.water` (lake depth) and `world.flow` (upstream rainfall).
    pub fn route(&self, world: &mut World, net: &mut FlowNetwork) {
        let p = &self.params;
        let n = world.len();
        net.filled.clear();
        net.filled.resize(n, f32::MAX);
        net.receiver.clear();
        net.receiver.resize(n, -1);
        net.order.clear();

        let mut heap: BinaryHeap<Entry> = BinaryHeap::with_capacity(n / 4);
        let mut queued = vec![false; n];
        // `filled` carries the epsilon tilt that guarantees drainage; that
        // tilt would otherwise show up as a millimetre of phantom water on
        // every flat plateau, so lake depth is measured against an untilted
        // copy of the same flood.
        let mut water_level = vec![f32::MAX; n];

        // Outlets: the sea, and the northern and southern walls (water that
        // reaches a wall leaves the world).
        for idx in 0..n {
            let (_, y) = world.coords(idx);
            let is_edge = y == 0 || y == world.height - 1;
            if world.elevation[idx] <= p.sea_level || is_edge {
                net.filled[idx] = world.elevation[idx];
                water_level[idx] = world.elevation[idx];
                queued[idx] = true;
                heap.push(Entry {
                    level: net.filled[idx],
                    idx: idx as u32,
                });
            }
        }

        while let Some(Entry { idx, .. }) = heap.pop() {
            let idx = idx as usize;
            net.order.push(idx as u32);
            let level = net.filled[idx];

            for (nidx, _, _) in world.neighbors8(idx) {
                if queued[nidx] {
                    continue;
                }
                queued[nidx] = true;
                // Either the cell's own height, or just above the water that
                // backed up to reach it.
                net.filled[nidx] = world.elevation[nidx].max(level + p.fill_epsilon);
                water_level[nidx] = world.elevation[nidx].max(water_level[idx]);
                heap.push(Entry {
                    level: net.filled[nidx],
                    idx: nidx as u32,
                });
            }
        }

        // Steepest descent on the filled surface. Because of the epsilon tilt
        // every non-outlet cell has somewhere to go.
        for idx in 0..n {
            if world.elevation[idx] <= p.sea_level {
                net.receiver[idx] = -1;
                world.water[idx] = 0.0;
                continue;
            }
            world.water[idx] = (water_level[idx] - world.elevation[idx]).max(0.0);

            let here = net.filled[idx];
            let mut best = -1i32;
            let mut best_drop = 0.0f32;
            for (nidx, dx, dy) in world.neighbors8(idx) {
                let dist = if dx != 0 && dy != 0 {
                    std::f32::consts::SQRT_2
                } else {
                    1.0
                };
                let drop = (here - net.filled[nidx]) / dist;
                if drop > best_drop {
                    best_drop = drop;
                    best = nidx as i32;
                }
            }
            net.receiver[idx] = best;
        }

        // Accumulate downstream. `order` ascends, so walking it backwards
        // always visits a cell before the cell it drains into.
        world.flow.copy_from_slice(&world.rainfall);
        for &idx in net.order.iter().rev() {
            let idx = idx as usize;
            let r = net.receiver[idx];
            if r >= 0 {
                world.flow[r as usize] += world.flow[idx];
            }
        }
    }
}

/// Rain by latitude: an equatorial peak, dry subtropics around 30 degrees, and
/// a secondary mid-latitude storm-track peak. Shape only - amplitude is tuned
/// by the params.
fn latitude_rainfall(lat: f32) -> f32 {
    let itcz = (-(lat / 11.0).powi(2)).exp();
    let storm_track = 0.55 * (-((lat.abs() - 52.0) / 16.0).powi(2)).exp();
    0.95 * itcz + storm_track
}

/// Prevailing surface wind: trade easterlies, mid-latitude westerlies, polar
/// easterlies. Only the x direction matters for the orographic term.
fn prevailing_wind(lat: f32) -> Vec2 {
    let a = lat.abs();
    if a < 30.0 {
        Vec2::new(-1.0, 0.0)
    } else if a < 60.0 {
        Vec2::new(1.0, 0.0)
    } else {
        Vec2::new(-1.0, 0.0)
    }
}

/// Min-heap entry: `BinaryHeap` is a max-heap, so the ordering is inverted.
#[derive(PartialEq)]
struct Entry {
    level: f32,
    idx: u32,
}

impl Eq for Entry {}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .level
            .total_cmp(&self.level)
            .then_with(|| other.idx.cmp(&self.idx))
    }
}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Distance between a cell and its receiver, for slope calculations.
pub fn step_distance(world: &World, from: usize, to: usize) -> f32 {
    let (fx, fy) = world.coords(from);
    let (tx, ty) = world.coords(to);
    let d = world.delta(tx as f32, ty as f32, fx as f32, fy as f32);
    let _ = OFFSETS8;
    d.length().max(1e-6)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_world(w: usize, h: usize, elev: f32) -> World {
        let mut world = World::new(w, h);
        world.elevation.fill(elev);
        world
    }

    fn routed(world: &mut World, hydro: &Hydrology) -> FlowNetwork {
        let mut net = FlowNetwork::new(world.len());
        hydro.rainfall(world);
        hydro.route(world, &mut net);
        net
    }

    fn hydro() -> Hydrology {
        Hydrology::new(HydrologyParams::default(), 7, 64)
    }

    #[test]
    fn filled_surface_never_dips_below_terrain() {
        let mut world = flat_world(64, 32, 1.0);
        for i in 0..world.len() {
            let (x, y) = world.coords(i);
            world.elevation[i] = 1.0 + ((x * 7 + y * 13) % 11) as f32 * 0.1;
        }
        let h = hydro();
        let net = routed(&mut world, &h);
        for i in 0..world.len() {
            assert!(
                net.filled[i] >= world.elevation[i] - 1e-6,
                "filled below terrain at {i}"
            );
        }
    }

    #[test]
    fn every_land_cell_drains_to_an_outlet() {
        let mut world = flat_world(64, 32, 2.0);
        for i in 0..world.len() {
            let (x, y) = world.coords(i);
            // Lumpy terrain with plenty of closed basins.
            world.elevation[i] =
                2.0 + ((x as f32 * 0.3).sin() * (y as f32 * 0.4).cos()) * 1.5;
        }
        let h = hydro();
        let net = routed(&mut world, &h);

        for start in 0..world.len() {
            let mut cur = start as i32;
            let mut hops = 0;
            while cur >= 0 {
                cur = net.receiver[cur as usize];
                hops += 1;
                assert!(hops <= world.len(), "cycle in the drainage network from {start}");
            }
        }
    }

    #[test]
    fn a_basin_fills_to_its_rim_and_holds_water() {
        // A plateau with a hole punched in the middle of it.
        let mut world = flat_world(48, 48, 3.0);
        for y in 20..28 {
            for x in 20..28 {
                let i = world.idx(x, y);
                world.elevation[i] = 0.5;
            }
        }
        let h = hydro();
        let net = routed(&mut world, &h);

        let centre = world.idx(24, 24);
        assert!(world.water[centre] > 2.0, "basin should hold water: {}", world.water[centre]);
        assert!(
            (net.filled[centre] - 3.0).abs() < 0.05,
            "basin should fill to the rim, got {}",
            net.filled[centre]
        );
        // Terrain outside the basin is not underwater.
        assert_eq!(world.water[world.idx(5, 5)], 0.0);
    }

    #[test]
    fn all_rainfall_reaches_a_sink() {
        let mut world = flat_world(48, 24, 1.0);
        for i in 0..world.len() {
            let (x, y) = world.coords(i);
            world.elevation[i] = 1.0 + (x as f32 * 0.2).sin() + (y as f32 * 0.3).cos();
        }
        let h = hydro();
        let net = routed(&mut world, &h);

        let total_rain: f64 = world.rainfall.iter().map(|r| *r as f64).sum();
        let at_sinks: f64 = (0..world.len())
            .filter(|i| net.receiver[*i] < 0)
            .map(|i| world.flow[i] as f64)
            .sum();
        assert!(
            (total_rain - at_sinks).abs() / total_rain < 1e-4,
            "water is not conserved: {total_rain} fell, {at_sinks} reached a sink"
        );
    }

    #[test]
    fn flow_grows_downstream() {
        let mut world = flat_world(64, 32, 1.0);
        for i in 0..world.len() {
            let (_, y) = world.coords(i);
            world.elevation[i] = 5.0 - y as f32 * 0.15; // slope toward the south wall
        }
        let h = hydro();
        let net = routed(&mut world, &h);
        for i in 0..world.len() {
            let r = net.receiver[i];
            if r >= 0 {
                assert!(
                    world.flow[r as usize] >= world.flow[i] - 1e-3,
                    "flow shrank downstream at {i}"
                );
            }
        }
    }

    #[test]
    fn drainage_crosses_the_seam() {
        // A valley floor running along x=0: everything drains over the seam.
        let mut world = flat_world(64, 32, 0.0);
        for i in 0..world.len() {
            let (x, y) = world.coords(i);
            let dx = crate::math::wrap_delta(x as f32, 0.0, 64.0).abs();
            world.elevation[i] = 1.0 + dx * 0.05 + (y as f32 - 16.0).abs() * 0.01;
        }
        let h = hydro();
        let net = routed(&mut world, &h);

        let west_of_seam = world.idx(63, 16);
        let mut cur = west_of_seam as i32;
        let mut crossed = false;
        let mut hops = 0;
        while cur >= 0 && hops < 1000 {
            let (x, _) = world.coords(cur as usize);
            if x <= 1 {
                crossed = true;
                break;
            }
            cur = net.receiver[cur as usize];
            hops += 1;
        }
        assert!(crossed, "water at x=63 never found the valley across the seam");
        assert!(net.receiver[west_of_seam] >= 0);
    }

    #[test]
    fn rainfall_is_wet_at_the_equator_and_dry_in_the_subtropics() {
        let mut world = flat_world(64, 180, -1.0); // all ocean: no orography
        let h = Hydrology::new(HydrologyParams::default(), 7, 64);
        h.rainfall(&mut world);

        let row_mean = |y: usize| {
            let mut s = 0.0;
            for x in 0..world.width {
                s += world.rainfall[world.idx(x, y)];
            }
            s / world.width as f32
        };
        let equator = row_mean(90);
        let subtropics = row_mean(120); // about 30 degrees south
        let pole = row_mean(178);
        assert!(equator > subtropics * 1.5, "equator {equator} vs subtropics {subtropics}");
        assert!(equator > pole, "equator {equator} vs pole {pole}");
    }

    #[test]
    fn mountains_cast_a_rain_shadow_on_their_downwind_side() {
        // A north/south ridge across every latitude band.
        let mut world = flat_world(64, 32, 0.2);
        for y in 0..world.height {
            for x in 28..32 {
                let i = world.idx(x, y);
                world.elevation[i] = 4.0;
            }
        }
        let h = hydro();
        h.rainfall(&mut world);
        let rain = |x: usize, y: usize| world.rainfall[world.idx(x, y)];

        // y=8 is about 42 degrees north: mid-latitude westerlies, so the
        // western flank climbs into the wind and the eastern flank is starved.
        assert!(
            rain(28, 8) > rain(32, 8) * 1.5,
            "westerlies: west flank {} should beat east flank {}",
            rain(28, 8),
            rain(32, 8)
        );

        // y=16 is near the equator: trade easterlies, so it is the other way
        // round. The wet flank follows the wind, not the map.
        assert!(
            rain(31, 16) > rain(27, 16) * 1.5,
            "trades: east flank {} should beat west flank {}",
            rain(31, 16),
            rain(27, 16)
        );
    }

    #[test]
    fn routing_is_deterministic() {
        let mut a = flat_world(48, 24, 1.0);
        for i in 0..a.len() {
            let (x, y) = a.coords(i);
            a.elevation[i] = (x as f32 * 0.11).sin() + (y as f32 * 0.17).cos();
        }
        let mut b = a.clone();
        let h = hydro();
        let na = routed(&mut a, &h);
        let nb = routed(&mut b, &h);
        assert_eq!(na.receiver, nb.receiver);
        assert_eq!(na.order, nb.order);
        assert_eq!(a.flow, b.flow);
    }
}
