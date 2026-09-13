//! The simulation loop (design doc section 18).
//!
//! One iteration is: move the plates, rebuild ownership, find and classify the
//! boundaries, deform the terrain through them, erupt volcanoes, then let
//! water rework the result. Every stage lives in its own module; this one only
//! decides the order and carries the buffers between them.

use crate::boundary::{self, Boundary, BoundaryKind, NearestBoundary};
use crate::erosion::{Erosion, ErosionBudget, ErosionParams};
use crate::hydrology::{FlowNetwork, Hydrology, HydrologyParams};
use crate::math::{lerp, smoothstep};
use crate::noise::CylinderNoise;
use crate::plate::{self, CrustField, Plate, PlateParams};
use crate::tectonics::{TectonicParams, Tectonics, mean_crust_age};
use crate::voronoi::{self, Warp};
use crate::world::{DEFAULT_HEIGHT, DEFAULT_WIDTH, World};

/// Clock based seed
fn clock_seed() -> [u64; 2] {
    let a = std::time::UNIX_EPOCH.elapsed().unwrap().as_millis();
    [(a >> 64) as u64, a as u64]
}

#[derive(Clone, Copy, Debug)]
pub struct SimulationParams {
    pub width: usize,
    pub height: usize,
    pub seed: u64,
    /// Geological time per iteration. Everything scales by this.
    pub dt: f32,
    /// Where the sea sits. This is the single source of truth: hydrology,
    /// erosion, the land/sea statistics and every renderer read it from here,
    /// so there is no way for two stages to disagree about what is underwater.
    ///
    /// Lowering it drains the world without touching the tectonics, which is
    /// the difference between a land-heavy world that still has trenches and
    /// island arcs, and one made of wall-to-wall continental crust.
    pub sea_level: f32,
    /// Normal-motion magnitude below which a boundary counts as transform.
    pub transform_threshold: f32,
    /// Run hydrology and erosion every N tectonic steps. Water works far
    /// faster than plates do, so it does not need to run every step.
    pub erosion_interval: u32,

    /// Starting elevations by crust type, before any tectonics.
    pub oceanic_elevation: f32,
    pub continental_elevation: f32,
    /// Amplitude and scale of the noise applied to the starting terrain.
    pub terrain_noise: f32,
    pub terrain_feature_size: f32,
    /// How far plate boundaries are bent away from straight Voronoi edges,
    /// in cells, and over what wavelength. Without this, every mountain belt
    /// is a straight line.
    pub boundary_warp: f32,
    pub boundary_warp_scale: f32,
    /// How much noise is mixed into the coastline. Zero gives smooth contours
    /// straight off the crust field; higher values give ragged coasts,
    /// peninsulas and offshore islands.
    pub coast_roughness: f32,
    /// Starting crust age spread, so ocean/ocean subduction has something to
    /// choose between on step one.
    pub initial_age_spread: f32,

    pub plates: PlateParams,
    pub tectonics: TectonicParams,
    pub hydrology: HydrologyParams,
    pub erosion: ErosionParams,
}

impl Default for SimulationParams {
    fn default() -> Self {
        SimulationParams {
            width: DEFAULT_WIDTH,
            height: DEFAULT_HEIGHT,
            seed: clock_seed()[1],
            dt: 1.0,
            sea_level: 0.0,
            transform_threshold: 0.04,
            erosion_interval: 2,

            oceanic_elevation: -3.6,
            continental_elevation: 1.15,
            terrain_noise: 1.75,
            terrain_feature_size: 55.0,
            coast_roughness: 0.42,
            boundary_warp: 22.0,
            boundary_warp_scale: 95.0,
            initial_age_spread: 60.0,

            plates: PlateParams::default(),
            tectonics: TectonicParams::default(),
            hydrology: HydrologyParams::default(),
            erosion: ErosionParams::default(),
        }
    }
}

/// What happened during one iteration - enough to tell whether a world went
/// wrong in the plates, the boundaries, the deformation or the water.
#[derive(Clone, Copy, Debug, Default)]
pub struct StepStats {
    pub step: u64,
    pub convergent: usize,
    pub divergent: usize,
    pub transform: usize,
    pub volcanoes: usize,
    pub min_elevation: f32,
    pub max_elevation: f32,
    pub land_fraction: f32,
    pub erosion: ErosionBudget,
}

pub struct Simulation {
    pub params: SimulationParams,
    pub world: World,
    pub plates: Vec<Plate>,
    pub boundaries: Vec<Boundary>,
    pub field: NearestBoundary,
    pub net: FlowNetwork,
    pub step: u64,
    pub stats: StepStats,
    /// Running erosion totals for the whole run. `stats.erosion` only covers
    /// the last step, and erosion does not run on every step.
    pub totals: ErosionBudget,

    tectonics: Tectonics,
    hydrology: Hydrology,
    erosion: Erosion,
    warp: Warp,
}

impl Simulation {
    pub fn new(params: SimulationParams) -> Self {
        let mut world = World::new(params.width, params.height);

        // Geography first, plates second: the crust field decides where
        // continents are, and the plates inherit their crust type from it.
        let crust = CrustField::new(&world, params.seed, params.plates.continental_fraction);
        let plates = plate::generate(&world, &params.plates, params.seed, &crust);
        let warp = Warp::new(
            &world,
            params.seed,
            params.boundary_warp,
            params.boundary_warp_scale,
        );
        voronoi::assign(&mut world, &plates, &warp);

        // Starting terrain. `shape` is the large-scale land/sea pattern with
        // enough noise mixed in that the coastline is ragged rather than a
        // clean contour; `relief` is the local texture that gives water
        // somewhere to run before erosion has carved anything.
        let shape = CylinderNoise::new(
            (params.seed as u32) ^ 0x_7E44,
            params.width,
            params.terrain_feature_size,
            5,
        );
        let relief = CylinderNoise::new(
            (params.seed as u32) ^ 0x_9C13,
            params.width,
            params.terrain_feature_size * 0.22,
            4,
        );

        // Starting crust age comes from its own noise field. Keying it to the
        // plate instead would stamp the Voronoi polygons into the sea floor
        // through age-driven subsidence.
        let age_field = CylinderNoise::new(
            (params.seed as u32) ^ 0x_A6E5,
            params.width,
            180.0,
            3,
        );

        for idx in 0..world.len() {
            let (x, y) = world.coords(idx);
            let (fx, fy) = (x as f32, y as f32);

            let coast = crust.continentalness(fx, fy) + params.coast_roughness * shape.get(fx, fy);
            let land = smoothstep(-0.04, 0.10, coast);
            let base = lerp(params.oceanic_elevation, params.continental_elevation, land);

            // Land is rougher than the sea floor, and the margin in between is
            // smoothed so shelves do not turn into cliffs.
            let texture = 0.6 * shape.get(fx, fy) + 0.4 * relief.get(fx, fy);
            world.elevation[idx] = base + texture * params.terrain_noise * (0.35 + 0.65 * land);
            world.crust_age[idx] = (0.5 + 0.5 * age_field.get(fx, fy)).clamp(0.0, 1.0)
                * params.initial_age_spread;
        }

        let tectonics = Tectonics::new(params.tectonics, params.seed, params.width);
        let hydrology =
            Hydrology::new(params.hydrology, params.sea_level, params.seed, params.width);
        let erosion = Erosion::new(params.erosion, params.sea_level, world.len());
        let boundaries = boundary::detect(&world, &plates, params.transform_threshold);
        let field = boundary::nearest_boundary(&world, &boundaries, tectonics.max_radius());
        let net = FlowNetwork::new(world.len());

        let mut sim = Simulation {
            params,
            world,
            plates,
            boundaries,
            field,
            net,
            step: 0,
            stats: StepStats::default(),
            totals: ErosionBudget::default(),
            tectonics,
            hydrology,
            erosion,
            warp,
        };
        // Give the initial world a drainage network so it can be rendered
        // before any steps have run.
        sim.water_cycle();
        sim.stats = sim.collect_stats(ErosionBudget::default());
        sim
    }

    /// One iteration of the sequence in design doc section 18.
    pub fn step(&mut self) {
        let dt = self.params.dt;

        // 1. Move the plates. 2. Rebuild ownership.
        plate::advance(&mut self.plates, &self.world, dt);
        voronoi::assign(&mut self.world, &self.plates, &self.warp);

        // 3-7. Find boundaries, take their normals, compare plate velocities,
        // classify, and work out how hard each one is pushing. All of that is
        // decided inside `detect`.
        self.boundaries =
            boundary::detect(&self.world, &self.plates, self.params.transform_threshold);
        self.field =
            boundary::nearest_boundary(&self.world, &self.boundaries, self.tectonics.max_radius());

        // 8-9, 12. Deform the terrain, and age the crust.
        let plate_age = mean_crust_age(&self.world, self.plates.len());
        self.tectonics.apply(
            &mut self.world,
            &self.plates,
            &self.boundaries,
            &self.field,
            &plate_age,
            dt,
        );

        // 10. Volcanism.
        let volcanoes = self.tectonics.apply_volcanism(
            &mut self.world,
            &self.plates,
            &self.boundaries,
            &plate_age,
            self.step,
            dt,
        );

        // 11. Water, on its own slower cadence.
        let mut budget = ErosionBudget::default();
        if self.params.erosion_interval > 0
            && self.step.is_multiple_of(self.params.erosion_interval as u64)
        {
            budget = self.water_cycle();
        }

        self.totals.eroded += budget.eroded;
        self.totals.deposited += budget.deposited;
        self.totals.exported += budget.exported;

        self.step += 1;
        self.stats = self.collect_stats(budget);
        self.stats.volcanoes = volcanoes;
    }

    pub fn run(&mut self, steps: u64) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// Rainfall, routing, erosion. Returns what the erosion moved.
    fn water_cycle(&mut self) -> ErosionBudget {
        self.hydrology.rainfall(&mut self.world);
        self.hydrology.route(&mut self.world, &mut self.net);
        self.erosion
            .apply(&mut self.world, &self.net, self.params.dt)
    }

    fn collect_stats(&self, erosion: ErosionBudget) -> StepStats {
        let mut stats = StepStats {
            step: self.step,
            erosion,
            min_elevation: f32::MAX,
            max_elevation: f32::MIN,
            ..Default::default()
        };
        for b in &self.boundaries {
            match b.kind {
                BoundaryKind::Convergent => stats.convergent += 1,
                BoundaryKind::Divergent => stats.divergent += 1,
                BoundaryKind::Transform => stats.transform += 1,
            }
        }
        let sea = self.params.sea_level;
        let mut land = 0usize;
        for &e in &self.world.elevation {
            stats.min_elevation = stats.min_elevation.min(e);
            stats.max_elevation = stats.max_elevation.max(e);
            if e > sea {
                land += 1;
            }
        }
        stats.land_fraction = land as f32 / self.world.len() as f32;
        stats
    }

    pub fn sea_level(&self) -> f32 {
        self.params.sea_level
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A reduced world for tests.
    ///
    /// Not tiny: the default deformation radii are absolute cell counts tuned
    /// for a 1024x512 world (design doc section 11), so in a very small world
    /// every rule reaches every cell and the plate interiors the model depends
    /// on do not exist. 256x128 keeps plates comfortably wider than the widest
    /// radius while still running in a fraction of a second.
    fn small(seed: u64) -> SimulationParams {
        SimulationParams {
            width: 256,
            height: 128,
            seed,
            plates: PlateParams {
                count: 12,
                ..Default::default()
            },
            ..Default::default()
        }
    }

    #[test]
    fn same_seed_gives_the_same_world() {
        let mut a = Simulation::new(small(4242));
        let mut b = Simulation::new(small(4242));
        a.run(25);
        b.run(25);
        assert_eq!(a.world.elevation, b.world.elevation);
        assert_eq!(a.world.plate_id, b.world.plate_id);
        assert_eq!(a.world.crust_age, b.world.crust_age);
        assert_eq!(a.world.flow, b.world.flow);
    }

    #[test]
    fn different_seeds_give_different_worlds() {
        let mut a = Simulation::new(small(1));
        let mut b = Simulation::new(small(2));
        a.run(10);
        b.run(10);
        assert_ne!(a.world.elevation, b.world.elevation);
    }

    #[test]
    fn a_long_run_stays_finite_and_bounded() {
        let mut sim = Simulation::new(small(9));
        sim.run(150);
        for (i, e) in sim.world.elevation.iter().enumerate() {
            assert!(e.is_finite(), "cell {i} went non-finite");
            assert!(
                (-30.0..30.0).contains(e),
                "cell {i} left any plausible range: {e}"
            );
        }
    }

    #[test]
    fn plates_keep_moving_and_ownership_keeps_changing() {
        let mut sim = Simulation::new(small(11));
        let start_positions: Vec<_> = sim.plates.iter().map(|p| p.position).collect();
        let start_owners = sim.world.plate_id.clone();
        sim.run(30);

        assert!(
            sim.plates
                .iter()
                .zip(&start_positions)
                .any(|(p, s)| p.position != *s),
            "plates never moved"
        );
        let changed = sim
            .world
            .plate_id
            .iter()
            .zip(&start_owners)
            .filter(|(a, b)| a != b)
            .count();
        assert!(changed > 0, "no cell ever changed hands");
    }

    #[test]
    fn all_three_boundary_types_show_up() {
        let mut sim = Simulation::new(small(3));
        sim.run(5);
        assert!(sim.stats.convergent > 0, "no convergent boundaries");
        assert!(sim.stats.divergent > 0, "no divergent boundaries");
        assert!(sim.stats.transform > 0, "no transform boundaries");
    }

    #[test]
    fn tectonics_build_relief_and_erosion_wears_it_down() {
        let mut sim = Simulation::new(small(17));
        let start_relief = sim.stats.max_elevation - sim.stats.min_elevation;
        sim.run(120);
        let relief = sim.stats.max_elevation - sim.stats.min_elevation;
        assert!(
            relief > start_relief,
            "120 steps of tectonics produced no new relief: {start_relief} -> {relief}"
        );
        assert!(sim.totals.eroded > 0.0, "erosion never ran");
        assert!(sim.totals.deposited > 0.0, "nothing was ever deposited");
    }

    #[test]
    fn sea_level_moves_the_coastline_without_touching_the_tectonics() {
        // Same seed and the same plates, but a lower sea: strictly more land,
        // and the terrain itself is not rebuilt around the new coastline.
        let mut normal = Simulation::new(small(7));
        let mut drained = Simulation::new(SimulationParams {
            sea_level: -1.5,
            ..small(7)
        });
        normal.run(40);
        drained.run(40);

        assert!(
            drained.stats.land_fraction > normal.stats.land_fraction,
            "lowering the sea should expose land: {} vs {}",
            drained.stats.land_fraction,
            normal.stats.land_fraction
        );
        assert_eq!(
            normal.world.plate_id, drained.world.plate_id,
            "sea level must not change which plate owns what"
        );
        assert_eq!(normal.sea_level(), 0.0);
        assert_eq!(drained.sea_level(), -1.5);
    }

    #[test]
    fn the_world_keeps_both_land_and_sea() {
        let mut sim = Simulation::new(small(23));
        sim.run(100);
        assert!(
            (0.05..0.95).contains(&sim.stats.land_fraction),
            "world is all one thing: land fraction {}",
            sim.stats.land_fraction
        );
    }

    #[test]
    fn crust_age_grows_but_spreading_centres_stay_young() {
        let mut sim = Simulation::new(small(31));
        sim.run(60);
        let mean: f32 =
            sim.world.crust_age.iter().sum::<f32>() / sim.world.len() as f32;
        let youngest = sim.world.crust_age.iter().cloned().fold(f32::MAX, f32::min);
        assert!(mean > 30.0, "crust should have aged: mean {mean}");
        assert!(youngest < 10.0, "nothing is being resurfaced: youngest {youngest}");
    }
}
