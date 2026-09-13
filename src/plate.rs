//! Tectonic plates: moving Voronoi sites (design doc section 3).

use crate::math::Vec2;
use crate::noise::{CylinderNoise, rng};
use crate::world::World;
use rand::RngExt as _;

/// Where continental crust sits, independent of where the plates are.
///
/// If crust type were read straight off plate ownership, every coastline would
/// be a Voronoi edge and the world would look like a polygon soup. Instead a
/// large-scale noise field decides what is continent and what is ocean; plates
/// then take their crust type from the field under their site, and the initial
/// terrain takes its shape from the same field. Continents and plates line up,
/// but coastlines are organic and a plate may straddle both (left open by
/// design doc section 26).
pub struct CrustField {
    noise: CylinderNoise,
    /// Field value above which crust is continental.
    pub threshold: f32,
}

impl CrustField {
    /// The threshold is picked so continental crust covers `continental_fraction`
    /// of the world, whatever the noise happens to look like for this seed.
    pub fn new(world: &World, seed: u64, continental_fraction: f32) -> Self {
        let noise = CylinderNoise::new((seed as u32) ^ 0x5EED, world.width, 260.0, 3);

        let mut samples: Vec<f32> = Vec::with_capacity(world.len() / 16 + 1);
        for y in (0..world.height).step_by(4) {
            for x in (0..world.width).step_by(4) {
                samples.push(noise.get(x as f32, y as f32));
            }
        }
        let fraction = continental_fraction.clamp(0.0, 1.0);
        let k = (((samples.len() - 1) as f32) * (1.0 - fraction)).round() as usize;
        let (_, nth, _) = samples.select_nth_unstable_by(k, f32::total_cmp);
        let threshold = *nth;

        CrustField { noise, threshold }
    }

    /// How continental a point is: positive on land-forming crust, negative on
    /// oceanic, and near zero along the continental margin.
    pub fn continentalness(&self, x: f32, y: f32) -> f32 {
        self.noise.get(x, y) - self.threshold
    }

    pub fn is_continental(&self, x: f32, y: f32) -> bool {
        self.continentalness(x, y) > 0.0
    }
}

pub type PlateId = u16;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CrustType {
    Oceanic,
    Continental,
}

impl CrustType {
    pub fn is_continental(self) -> bool {
        self == CrustType::Continental
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Plate {
    pub id: PlateId,
    /// Voronoi site, in cell coordinates.
    pub position: Vec2,
    /// Cells per simulation step (before `dt`).
    pub velocity: Vec2,
    pub crust_type: CrustType,
}

#[derive(Clone, Copy, Debug)]
pub struct PlateParams {
    pub count: usize,
    /// Fraction of plates that are continental.
    pub continental_fraction: f32,
    /// Plate speed range, in cells per unit time.
    pub speed_min: f32,
    pub speed_max: f32,
    /// Candidates per site for Mitchell's best-candidate spacing.
    pub spacing_candidates: usize,
}

impl Default for PlateParams {
    fn default() -> Self {
        PlateParams {
            count: 32,
            continental_fraction: 0.4,
            speed_min: 0.05,
            speed_max: 0.35,
            spacing_candidates: 24,
        }
    }
}

/// Place plate sites, give them velocities, and read off their crust types.
///
/// Sites are spread with Mitchell's best-candidate sampling so plates come out
/// roughly evenly sized instead of clumping. Crust type comes from the
/// [`CrustField`] under each site, which is what keeps plates consistent with
/// the geography the world was built from.
pub fn generate(world: &World, params: &PlateParams, seed: u64, crust: &CrustField) -> Vec<Plate> {
    let mut r = rng(seed, 0x91A7);
    let w = world.width as f32;
    let h = world.height as f32;

    let mut positions: Vec<Vec2> = Vec::with_capacity(params.count);
    for _ in 0..params.count {
        let mut best = Vec2::ZERO;
        let mut best_dist = -1.0f32;
        let candidates = params.spacing_candidates.max(1);
        for _ in 0..candidates {
            let c = Vec2::new(r.random_range(0.0..w), r.random_range(0.0..h));
            let d = positions
                .iter()
                .map(|p| world.dist2(c.x, c.y, p.x, p.y))
                .fold(f32::MAX, f32::min);
            if positions.is_empty() || d > best_dist {
                best_dist = d;
                best = c;
            }
        }
        positions.push(best);
    }

    positions
        .into_iter()
        .enumerate()
        .map(|(i, position)| {
            let angle = r.random_range(0.0..std::f32::consts::TAU);
            let speed =
                r.random_range(params.speed_min..params.speed_max.max(params.speed_min + 1e-6));
            Plate {
                id: i as PlateId,
                position,
                velocity: Vec2::new(angle.cos(), angle.sin()) * speed,
                crust_type: if crust.is_continental(position.x, position.y) {
                    CrustType::Continental
                } else {
                    CrustType::Oceanic
                },
            }
        })
        .collect()
}

/// Advance plate sites by one step.
///
/// East/west wraps. The N/S walls reflect: a site that would walk off the top
/// or bottom of the world bounces, which keeps plates in play instead of
/// letting them pile up against the wall (design doc section 26 left this
/// open; reflection is the least surprising choice).
pub fn advance(plates: &mut [Plate], world: &World, dt: f32) {
    let h = world.height as f32 - 1.0;
    for plate in plates.iter_mut() {
        let mut next = plate.position + plate.velocity * dt;
        if next.y < 0.0 {
            next.y = -next.y;
            plate.velocity.y = -plate.velocity.y;
        } else if next.y > h {
            next.y = 2.0 * h - next.y;
            plate.velocity.y = -plate.velocity.y;
        }
        plate.position = world.contain(next);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn world() -> World {
        World::new(256, 128)
    }

    fn field(world: &World, seed: u64, fraction: f32) -> CrustField {
        CrustField::new(world, seed, fraction)
    }

    fn make_plates(world: &World, params: &PlateParams, seed: u64) -> Vec<Plate> {
        let f = field(world, seed, params.continental_fraction);
        generate(world, params, seed, &f)
    }

    #[test]
    fn generation_is_deterministic() {
        let w = world();
        let p = PlateParams::default();
        let a = make_plates(&w, &p, 1234);
        let b = make_plates(&w, &p, 1234);
        let c = make_plates(&w, &p, 5678);
        assert_eq!(a.len(), p.count);
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.position, y.position);
            assert_eq!(x.velocity, y.velocity);
            assert_eq!(x.crust_type, y.crust_type);
        }
        assert!(a.iter().zip(c.iter()).any(|(x, y)| x.position != y.position));
    }

    #[test]
    fn the_crust_field_covers_the_requested_fraction_of_the_world() {
        let w = world();
        for fraction in [0.2f32, 0.4, 0.65] {
            let f = field(&w, 9, fraction);
            let mut land = 0;
            for y in 0..w.height {
                for x in 0..w.width {
                    if f.is_continental(x as f32, y as f32) {
                        land += 1;
                    }
                }
            }
            let actual = land as f32 / w.len() as f32;
            assert!(
                (actual - fraction).abs() < 0.03,
                "asked for {fraction} continental, got {actual}"
            );
        }
    }

    #[test]
    fn plates_take_their_crust_from_the_field_under_them() {
        let w = world();
        let params = PlateParams {
            count: 40,
            continental_fraction: 0.35,
            ..PlateParams::default()
        };
        let f = field(&w, 9, params.continental_fraction);
        let plates = generate(&w, &params, 9, &f);
        for p in &plates {
            assert_eq!(
                p.crust_type.is_continental(),
                f.is_continental(p.position.x, p.position.y),
                "plate {} disagrees with the crust under it",
                p.id
            );
        }
        // Roughly - not exactly - the requested share, since plates are a
        // sample of the field rather than a partition of it.
        let cont = plates.iter().filter(|p| p.crust_type.is_continental()).count();
        assert!((6..=22).contains(&cont), "{cont} of 40 plates continental");
    }

    #[test]
    fn sites_land_inside_the_world() {
        let w = world();
        for p in make_plates(&w, &PlateParams::default(), 3) {
            assert!(p.position.x >= 0.0 && p.position.x < w.width as f32);
            assert!(p.position.y >= 0.0 && p.position.y <= w.height as f32 - 1.0);
        }
    }

    #[test]
    fn movement_wraps_east_west_and_bounces_off_walls() {
        let w = world();
        let mut plates = vec![
            Plate {
                id: 0,
                position: Vec2::new(255.5, 10.0),
                velocity: Vec2::new(1.0, 0.0),
                crust_type: CrustType::Oceanic,
            },
            Plate {
                id: 1,
                position: Vec2::new(10.0, 1.0),
                velocity: Vec2::new(0.0, -4.0),
                crust_type: CrustType::Oceanic,
            },
        ];
        advance(&mut plates, &w, 1.0);
        assert!(plates[0].position.x < 1.0, "should have wrapped: {:?}", plates[0].position);
        assert!(plates[1].position.y >= 0.0);
        assert!(plates[1].velocity.y > 0.0, "wall should reflect velocity");
    }

    #[test]
    fn best_candidate_spacing_beats_pure_random() {
        let w = world();
        let spaced = make_plates(&w, &PlateParams { spacing_candidates: 32, ..Default::default() }, 11);
        let clumped = make_plates(&w, &PlateParams { spacing_candidates: 1, ..Default::default() }, 11);
        let min_gap = |ps: &Vec<Plate>| {
            let mut m = f32::MAX;
            for i in 0..ps.len() {
                for j in (i + 1)..ps.len() {
                    m = m.min(w.dist2(ps[i].position.x, ps[i].position.y, ps[j].position.x, ps[j].position.y));
                }
            }
            m
        };
        assert!(min_gap(&spaced) > min_gap(&clumped));
    }
}
