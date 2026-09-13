//! Reading a finished world.
//!
//! The simulation's fields are all public, and bulk consumers (building a tile
//! map, exporting a heightmap) should iterate them directly - that is the fast
//! path. This module is for the other case: asking everything about one place,
//! which otherwise means indexing five parallel arrays and joining plate
//! ownership against the plate list by hand.
//!
//! Nothing here mutates or advances anything, so a [`Simulation`] that has
//! finished running is a perfectly good read-only world model to keep around
//! for the life of a game.

use crate::boundary::{BoundaryKind, NO_SOURCE};
use crate::math::{Vec2, percentile};
use crate::plate::{CrustType, PlateId};
use crate::simulation::Simulation;

/// Everything the simulation knows about one cell.
///
/// Units are the simulation's own and are not metres. Sea level is 0.0; ocean
/// basins bottom out around -7 and the highest peaks reach about +8, so one
/// unit is roughly a kilometre if you want a mental scale.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sample {
    pub x: usize,
    pub y: usize,

    /// Height relative to sea level.
    pub elevation: f32,
    /// The world's sea level, repeated here so a `Sample` is self-contained.
    pub sea_level: f32,
    /// Latitude in degrees: +90 at the north wall, -90 at the south.
    pub latitude: f32,

    /// Crust type of the plate that currently holds this cell. This is what
    /// the tectonic rules key off, and it changes as plates move over the cell.
    pub crust: CrustType,
    /// Which plate holds the cell, and where that plate is heading, in cells
    /// per unit of simulation time.
    pub plate: PlateId,
    pub plate_velocity: Vec2,
    /// Time since this crust was created. Crust made at a spreading centre
    /// starts at zero, so this rises with distance from a mid-ocean ridge.
    pub crust_age: f32,

    /// Rain falling on this cell. Relative, not millimetres: compare cells to
    /// each other rather than reading an absolute figure off it.
    pub rainfall: f32,
    /// Depth of standing water above the terrain - a lake, if positive.
    pub lake_depth: f32,
    /// Upstream rainfall arriving here. Large values are rivers; see
    /// [`Simulation::river_threshold`].
    pub flow: f32,
    /// Sediment deposited here over the world's history.
    pub sediment: f32,

    /// The nearest plate boundary, if one is close enough to be shaping this
    /// cell at all. `None` means the cell is plate interior.
    pub boundary: Option<BoundaryProximity>,
}

/// The plate boundary currently shaping a cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryProximity {
    pub kind: BoundaryKind,
    /// Distance in cells to the boundary.
    pub distance: f32,
    /// How hard this boundary is working, from the relative plate motion.
    pub strength: f32,
    /// The two plates that meet here.
    pub plates: (PlateId, PlateId),
}

impl Sample {
    pub fn is_land(&self) -> bool {
        self.elevation > self.sea_level
    }

    pub fn is_sea(&self) -> bool {
        !self.is_land()
    }

    /// Depth below sea level, or 0 on land.
    pub fn depth(&self) -> f32 {
        (self.sea_level - self.elevation).max(0.0)
    }

    /// Height above sea level, or 0 at sea.
    pub fn altitude(&self) -> f32 {
        (self.elevation - self.sea_level).max(0.0)
    }

    pub fn is_lake(&self) -> bool {
        self.is_land() && self.lake_depth > 0.02
    }

    pub fn is_continental(&self) -> bool {
        self.crust.is_continental()
    }
}

impl Simulation {
    /// Everything known about one cell.
    ///
    /// Panics if the coordinates are outside the world; use
    /// [`Simulation::sample_at`] for coordinates that may need wrapping.
    pub fn sample(&self, x: usize, y: usize) -> Sample {
        let w = &self.world;
        let idx = w.idx(x, y);
        let plate = w.plate_id[idx];
        let p = &self.plates[plate as usize];

        let boundary = match self.field.source[idx] {
            NO_SOURCE => None,
            src => {
                let b = &self.boundaries[src as usize];
                Some(BoundaryProximity {
                    kind: b.kind,
                    distance: self.field.dist[idx],
                    strength: b.strength,
                    plates: (b.plate_a, b.plate_b),
                })
            }
        };

        Sample {
            x,
            y,
            elevation: w.elevation[idx],
            sea_level: self.sea_level(),
            latitude: w.latitude(y),
            crust: p.crust_type,
            plate,
            plate_velocity: p.velocity,
            crust_age: w.crust_age[idx],
            rainfall: w.rainfall[idx],
            lake_depth: w.water[idx],
            flow: w.flow[idx],
            sediment: w.sediment[idx],
            boundary,
        }
    }

    /// Like [`Simulation::sample`], but takes signed coordinates: x wraps
    /// around the cylinder, and a y off the north or south wall gives `None`.
    ///
    /// This is the one to use when a game walks off the edge of a region.
    pub fn sample_at(&self, x: i64, y: i64) -> Option<Sample> {
        if y < 0 || y >= self.world.height as i64 {
            return None;
        }
        Some(self.sample(self.world.wrap_x(x), y as usize))
    }

    /// Elevation at a continuous position, bilinearly interpolated.
    ///
    /// For a game whose coordinates are finer than one cell. Only elevation is
    /// interpolated: crust type and plate id are categorical, and interpolating
    /// rainfall or flow would invent water that is not there.
    pub fn elevation_at(&self, x: f32, y: f32) -> f32 {
        let w = &self.world;
        let y = y.clamp(0.0, w.height as f32 - 1.0);
        let (x0, y0) = (x.floor(), y.floor());
        let (tx, ty) = (x - x0, y - y0);

        let xa = w.wrap_x(x0 as i64);
        let xb = w.wrap_x(x0 as i64 + 1);
        let ya = w.clamp_y(y0 as i64);
        let yb = w.clamp_y(y0 as i64 + 1);

        let e = |cx: usize, cy: usize| w.elevation[w.idx(cx, cy)];
        let top = e(xa, ya) + (e(xb, ya) - e(xa, ya)) * tx;
        let bottom = e(xa, yb) + (e(xb, yb) - e(xa, yb)) * tx;
        top + (bottom - top) * ty
    }

    /// Is this cell above sea level?
    pub fn is_land(&self, x: usize, y: usize) -> bool {
        self.world.elevation[self.world.idx(x, y)] > self.sea_level()
    }

    /// Flow accumulation above which a cell is worth calling a river.
    ///
    /// Derived from this world's own distribution rather than being a fixed
    /// number, since total flow scales with world size and rainfall. The top
    /// half percent of cells by flow, which is about what reads as a river
    /// network on the map.
    ///
    /// O(cells), so hoist it out of a loop rather than calling it per cell.
    pub fn river_threshold(&self) -> f32 {
        percentile(&self.world.flow, 0.995).max(1.0)
    }

    /// Cells whose flow is at or above `threshold` and that are above sea
    /// level - the river network, as cell indices.
    pub fn rivers(&self, threshold: f32) -> impl Iterator<Item = Sample> + '_ {
        let sea = self.sea_level();
        (0..self.world.len()).filter_map(move |idx| {
            if self.world.flow[idx] >= threshold && self.world.elevation[idx] > sea {
                let (x, y) = self.world.coords(idx);
                Some(self.sample(x, y))
            } else {
                None
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plate::PlateParams;
    use crate::simulation::SimulationParams;

    fn world() -> Simulation {
        let mut sim = Simulation::new(SimulationParams {
            width: 256,
            height: 128,
            seed: 88,
            plates: PlateParams {
                count: 12,
                ..Default::default()
            },
            ..Default::default()
        });
        sim.run(60);
        sim
    }

    #[test]
    fn a_sample_agrees_with_the_raw_fields() {
        let sim = world();
        for (x, y) in [(0usize, 0usize), (17, 43), (255, 127), (128, 64)] {
            let s = sim.sample(x, y);
            let idx = sim.world.idx(x, y);
            assert_eq!(s.elevation, sim.world.elevation[idx]);
            assert_eq!(s.plate, sim.world.plate_id[idx]);
            assert_eq!(s.crust, sim.plates[s.plate as usize].crust_type);
            assert_eq!(s.crust_age, sim.world.crust_age[idx]);
            assert_eq!(s.rainfall, sim.world.rainfall[idx]);
            assert_eq!(s.flow, sim.world.flow[idx]);
            assert_eq!(s.lake_depth, sim.world.water[idx]);
            assert_eq!(s.latitude, sim.world.latitude(y));
        }
    }

    #[test]
    fn land_and_sea_predicates_are_consistent() {
        let sim = world();
        for idx in (0..sim.world.len()).step_by(97) {
            let (x, y) = sim.world.coords(idx);
            let s = sim.sample(x, y);
            assert_eq!(s.is_land(), !s.is_sea());
            assert_eq!(s.is_land(), sim.is_land(x, y));
            if s.is_land() {
                assert!(s.depth() == 0.0 && s.altitude() > 0.0);
            } else {
                assert!(s.altitude() == 0.0 && s.depth() >= 0.0);
            }
        }
    }

    #[test]
    fn sampling_wraps_east_west_and_stops_at_the_walls() {
        let sim = world();
        let w = sim.world.width as i64;
        assert_eq!(sim.sample_at(-1, 10), Some(sim.sample(255, 10)));
        assert_eq!(sim.sample_at(w, 10), Some(sim.sample(0, 10)));
        assert_eq!(sim.sample_at(w + 5, 10), Some(sim.sample(5, 10)));
        assert_eq!(sim.sample_at(10, -1), None);
        assert_eq!(sim.sample_at(10, 128), None);
    }

    #[test]
    fn interpolated_elevation_matches_cells_and_blends_between_them() {
        let sim = world();
        for (x, y) in [(3usize, 9usize), (100, 50)] {
            assert!(
                (sim.elevation_at(x as f32, y as f32) - sim.sample(x, y).elevation).abs() < 1e-4,
                "interpolation should be exact on cell centres"
            );
        }

        // Halfway between two cells is halfway between their heights.
        let a = sim.sample(40, 40).elevation;
        let b = sim.sample(41, 40).elevation;
        let mid = sim.elevation_at(40.5, 40.0);
        assert!((mid - (a + b) * 0.5).abs() < 1e-4, "{mid} is not between {a} and {b}");
    }

    #[test]
    fn interpolation_wraps_across_the_seam() {
        let sim = world();
        let west = sim.sample(255, 60).elevation;
        let east = sim.sample(0, 60).elevation;
        let mid = sim.elevation_at(255.5, 60.0);
        assert!((mid - (west + east) * 0.5).abs() < 1e-4);
    }

    #[test]
    fn boundary_proximity_is_reported_only_near_boundaries() {
        let sim = world();
        let mut near = 0;
        let mut interior = 0;
        for idx in 0..sim.world.len() {
            let (x, y) = sim.world.coords(idx);
            match sim.sample(x, y).boundary {
                Some(b) => {
                    near += 1;
                    assert!(b.distance <= sim.field.max_radius);
                    assert!(b.strength >= 0.0);
                    assert_ne!(b.plates.0, b.plates.1);
                }
                None => interior += 1,
            }
        }
        assert!(near > 0, "no cell is near a boundary");
        assert!(interior > 0, "no cell is plate interior");
    }

    #[test]
    fn rivers_are_land_cells_above_the_threshold() {
        let sim = world();
        let threshold = sim.river_threshold();
        let rivers: Vec<Sample> = sim.rivers(threshold).collect();
        assert!(!rivers.is_empty(), "world has no rivers at all");
        for r in &rivers {
            assert!(r.is_land());
            assert!(r.flow >= threshold);
        }
        // A threshold in the top half percent should select a small minority.
        assert!(rivers.len() < sim.world.len() / 20);
    }

    /// A finished world is meant to be handed to a game and kept, possibly
    /// across threads, so this needs to hold.
    #[test]
    fn a_simulation_can_be_shared_between_threads() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<Simulation>();
        assert_send_sync::<Sample>();
    }
}
