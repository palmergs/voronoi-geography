//! Plate ownership: which plate's Voronoi cell each terrain point falls in.
//!
//! Recomputed every step, since the sites move (design doc section 5). The
//! distance metric is cylindrical, so the diagram wraps across the x seam.
//!
//! The brute-force form is O(cells x plates). At 1024x512 with ~32 plates that
//! is ~17M distance tests, which `rayon` chews through in a few milliseconds -
//! not worth a fancier structure until the plate count grows a lot.

use crate::noise::CylinderNoise;
use crate::plate::{Plate, PlateId};
use crate::world::World;
use rayon::prelude::*;

/// A fixed distortion of the plane, applied before the Voronoi lookup.
///
/// Without it, plate boundaries are exactly the perpendicular bisectors of the
/// sites, so every mountain belt comes out as a dead straight line and the
/// polygon pattern shows through the finished world. Warping the lookup
/// position by smooth noise bends the boundaries into something that looks
/// like a rift or a suture, at no cost to the rest of the model: the diagram
/// is still deterministic, still tiles across the seam, and boundary normals
/// are still the site-to-site direction to within the warp's own gradient.
///
/// The offsets are precomputed once, since the distortion belongs to the
/// world rather than to any particular step.
pub struct Warp {
    dx: Vec<f32>,
    dy: Vec<f32>,
}

impl Warp {
    /// `amplitude` is in cells; `feature_size` is how wide one wobble is.
    pub fn new(world: &World, seed: u64, amplitude: f32, feature_size: f32) -> Self {
        let nx = CylinderNoise::new((seed as u32) ^ 0x_1A2B, world.width, feature_size, 3);
        let ny = CylinderNoise::new((seed as u32) ^ 0x_3C4D, world.width, feature_size, 3);
        let mut dx = vec![0.0; world.len()];
        let mut dy = vec![0.0; world.len()];
        for idx in 0..world.len() {
            let (x, y) = world.coords(idx);
            dx[idx] = nx.get(x as f32, y as f32) * amplitude;
            dy[idx] = ny.get(x as f32, y as f32) * amplitude;
        }
        Warp { dx, dy }
    }

    /// A warp that does nothing, for tests and for the unwarped case.
    pub fn none(world: &World) -> Self {
        Warp {
            dx: vec![0.0; world.len()],
            dy: vec![0.0; world.len()],
        }
    }
}

/// Fill `world.plate_id` from the current plate sites.
pub fn assign(world: &mut World, plates: &[Plate], warp: &Warp) {
    assert!(!plates.is_empty(), "cannot build a Voronoi diagram with no sites");

    let width = world.width;
    let w = width as f32;
    let half_w = w * 0.5;

    world
        .plate_id
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            let row_y = y as f32;
            for (x, cell) in row.iter_mut().enumerate() {
                let idx = y * width + x;
                let fx = x as f32 + warp.dx[idx];
                let fy = row_y + warp.dy[idx];
                let mut best = f32::MAX;
                let mut best_id: PlateId = 0;
                for plate in plates.iter() {
                    let mut dx = (fx - plate.position.x).abs();
                    if dx > half_w {
                        dx = w - dx;
                    }
                    let dy = fy - plate.position.y;
                    let d2 = dx * dx + dy * dy;
                    if d2 < best {
                        best = d2;
                        best_id = plate.id;
                    }
                }
                *cell = best_id;
            }
        });
}

#[cfg(test)]
mod tests {
    // Note: the nearest-site tests deliberately use `Warp::none`, since a warp
    // is exactly a controlled violation of the nearest-site property.

    use super::*;
    use crate::math::Vec2;
    use crate::plate::{CrustField, CrustType, PlateParams, generate};

    fn plate_at(id: PlateId, x: f32, y: f32) -> Plate {
        Plate {
            id,
            position: Vec2::new(x, y),
            velocity: Vec2::ZERO,
            crust_type: CrustType::Oceanic,
        }
    }

    #[test]
    fn every_cell_goes_to_its_nearest_site() {
        let mut world = World::new(64, 32);
        let plates = generate(&world, &PlateParams { count: 8, ..Default::default() }, 5, &CrustField::new(&world, 5, 0.4));
        let flat = Warp::none(&world);
        assign(&mut world, &plates, &flat);

        for y in 0..world.height {
            for x in 0..world.width {
                let owner = world.plate_id[world.idx(x, y)] as usize;
                let owner_d = world.dist2(
                    x as f32,
                    y as f32,
                    plates[owner].position.x,
                    plates[owner].position.y,
                );
                for p in &plates {
                    let d = world.dist2(x as f32, y as f32, p.position.x, p.position.y);
                    assert!(owner_d <= d + 1e-4, "cell ({x},{y}) claimed by the wrong plate");
                }
            }
        }
    }

    #[test]
    fn ownership_wraps_across_the_seam() {
        let mut world = World::new(64, 8);
        // One site just east of the seam, one on the far side of the map.
        let plates = vec![plate_at(0, 1.0, 4.0), plate_at(1, 32.0, 4.0)];
        let flat = Warp::none(&world);
        assign(&mut world, &plates, &flat);
        // x=63 is 2 cells from site 0 the short way, 31 from site 1.
        assert_eq!(world.plate_id[world.idx(63, 4)], 0);
        assert_eq!(world.plate_id[world.idx(40, 4)], 1);
    }

    #[test]
    fn sites_own_their_own_cell() {
        let mut world = World::new(128, 64);
        let plates = generate(&world, &PlateParams { count: 12, ..Default::default() }, 77, &CrustField::new(&world, 77, 0.4));
        let flat = Warp::none(&world);
        assign(&mut world, &plates, &flat);
        for p in &plates {
            let idx = world.idx(p.position.x.round() as usize % 128, p.position.y.round() as usize);
            assert_eq!(world.plate_id[idx], p.id);
        }
    }

    #[test]
    fn a_warp_bends_boundaries_without_breaking_the_diagram() {
        let mut plain = World::new(128, 64);
        let mut warped = World::new(128, 64);
        let crust = CrustField::new(&plain, 12, 0.4);
        let plates = generate(&plain, &PlateParams { count: 10, ..Default::default() }, 12, &crust);

        let flat = Warp::none(&plain);

        assign(&mut plain, &plates, &flat);
        let bent = Warp::new(&warped, 12, 9.0, 40.0);
        assign(&mut warped, &plates, &bent);

        // Ownership changes, but every cell still belongs to a real plate and
        // every plate still holds territory.
        assert_ne!(plain.plate_id, warped.plate_id);
        for id in &warped.plate_id {
            assert!((*id as usize) < plates.len());
        }
        for p in &plates {
            assert!(
                warped.plate_id.contains(&p.id),
                "plate {} lost all of its territory to the warp",
                p.id
            );
        }
    }

    #[test]
    fn warped_boundaries_are_longer_than_straight_ones() {
        // A wigglier boundary means more boundary cells: that is the whole
        // point of the warp.
        let mut plain = World::new(256, 128);
        let mut warped = World::new(256, 128);
        let crust = CrustField::new(&plain, 4, 0.4);
        let plates = generate(&plain, &PlateParams { count: 14, ..Default::default() }, 4, &crust);
        let flat = Warp::none(&plain);
        assign(&mut plain, &plates, &flat);
        let bent = Warp::new(&warped, 4, 10.0, 45.0);
        assign(&mut warped, &plates, &bent);

        let edge_cells = |w: &World| {
            (0..w.len())
                .filter(|i| {
                    w.neighbors8(*i)
                        .any(|(n, _, _)| w.plate_id[n] != w.plate_id[*i])
                })
                .count()
        };
        assert!(
            edge_cells(&warped) > edge_cells(&plain),
            "the warp should lengthen boundaries, not shorten them"
        );
    }

    #[test]
    fn the_warp_is_deterministic() {
        let world = World::new(64, 32);
        let a = Warp::new(&world, 77, 8.0, 30.0);
        let b = Warp::new(&world, 77, 8.0, 30.0);
        assert_eq!(a.dx, b.dx);
        assert_eq!(a.dy, b.dy);
    }
}
