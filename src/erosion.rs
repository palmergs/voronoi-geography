//! Erosion and deposition (design doc section 17).
//!
//! Three processes, run in order on the drainage network hydrology just built:
//!
//! 1. **Stream-power incision** - rivers cut down in proportion to discharge
//!    and slope. This is what carves valleys into the tectonic uplift, and
//!    what makes big rivers cut deeper than small ones.
//! 2. **Sediment transport** - what a river cuts, it carries downstream and
//!    drops progressively, dumping the rest into lakes and the sea. Deltas and
//!    infilled basins come from here.
//! 3. **Hillslope diffusion** - slow creep that rounds off ridges and stops
//!    the heightfield from developing single-cell spikes.
//!
//! Incision is capped at a fraction of the drop to the downstream cell, which
//! is what keeps the explicit time stepping stable no matter how large the
//! uplift or discharge gets.

use crate::hydrology::{FlowNetwork, step_distance};
use crate::world::World;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
pub struct ErosionParams {
    /// Stream-power coefficient.
    pub incision: f32,
    /// Discharge exponent (m) and slope exponent (n) in `k * A^m * S^n`.
    pub area_exponent: f32,
    pub slope_exponent: f32,
    /// Stability clamp: never cut more than this share of the drop downstream.
    pub max_incision_fraction: f32,
    /// Share of the carried sediment load dropped per cell on land.
    pub deposition: f32,
    /// Hillslope creep coefficient.
    pub diffusion: f32,
    /// Creep is slower under water, where there is no rain to drive it.
    pub submarine_diffusion_scale: f32,
}

impl Default for ErosionParams {
    fn default() -> Self {
        ErosionParams {
            incision: 0.0016,
            area_exponent: 0.5,
            slope_exponent: 1.0,
            max_incision_fraction: 0.3,
            deposition: 0.08,
            diffusion: 0.022,
            submarine_diffusion_scale: 0.25,
        }
    }
}

/// What one erosion step moved, for reporting and for tests.
#[derive(Clone, Copy, Debug, Default)]
pub struct ErosionBudget {
    pub eroded: f64,
    pub deposited: f64,
    /// Sediment carried off the map at a wall, or lost to the deep ocean.
    pub exported: f64,
}

pub struct Erosion {
    pub params: ErosionParams,
    /// Owned by the simulation, not by this stage - see `SimulationParams`.
    sea_level: f32,
    /// Scratch, kept across steps to avoid reallocating every iteration.
    cut: Vec<f32>,
    carry: Vec<f32>,
    scratch: Vec<f32>,
}

impl Erosion {
    pub fn new(params: ErosionParams, sea_level: f32, n: usize) -> Self {
        Erosion {
            params,
            sea_level,
            cut: vec![0.0; n],
            carry: vec![0.0; n],
            scratch: vec![0.0; n],
        }
    }

    pub fn apply(&mut self, world: &mut World, net: &FlowNetwork, dt: f32) -> ErosionBudget {
        let mut budget = ErosionBudget::default();
        self.incise(world, net, dt, &mut budget);
        self.transport(world, net, dt, &mut budget);
        self.diffuse(world, dt);
        budget
    }

    /// Stream-power incision: `k * discharge^m * slope^n`, clamped for stability.
    fn incise(&mut self, world: &mut World, net: &FlowNetwork, dt: f32, budget: &mut ErosionBudget) {
        let p = &self.params;
        let sea = self.sea_level;
        let elevation = &world.elevation;
        let flow = &world.flow;
        let water = &world.water;

        self.cut
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, cut)| {
                *cut = 0.0;
                let r = net.receiver[idx];
                // No rivers below sea level, and none inside a lake - still
                // water does not incise.
                if r < 0 || elevation[idx] <= sea || water[idx] > 0.0 {
                    return;
                }
                let r = r as usize;
                let drop = elevation[idx] - elevation[r];
                if drop <= 0.0 {
                    return;
                }
                let slope = drop / step_distance(world, idx, r);
                let e = p.incision
                    * flow[idx].max(0.0).powf(p.area_exponent)
                    * slope.powf(p.slope_exponent)
                    * dt;
                *cut = e.min(drop * p.max_incision_fraction);
            });

        for idx in 0..world.len() {
            let c = self.cut[idx];
            if c > 0.0 {
                world.elevation[idx] -= c;
                budget.eroded += c as f64;
            }
        }
    }

    /// Carry the cut material downstream, dropping a share of it on the way.
    ///
    /// `order` ascends by filled elevation, so walking it backwards visits
    /// every cell before the cell it drains into - one pass, no iteration.
    fn transport(
        &mut self,
        world: &mut World,
        net: &FlowNetwork,
        dt: f32,
        budget: &mut ErosionBudget,
    ) {
        let p = &self.params;
        self.carry.copy_from_slice(&self.cut);

        for &idx in net.order.iter().rev() {
            let idx = idx as usize;
            let load = self.carry[idx];
            if load <= 0.0 {
                continue;
            }
            let r = net.receiver[idx];

            // A sink swallows everything: lakes silt up, river mouths build
            // out into the sea.
            if r < 0 {
                world.elevation[idx] += load;
                world.sediment[idx] += load;
                budget.deposited += load as f64;
                self.carry[idx] = 0.0;
                continue;
            }

            // Slack water (a lake) drops its whole load; moving water drops a
            // share of it.
            let share = if world.water[idx] > 0.0 {
                1.0
            } else {
                (p.deposition * dt).clamp(0.0, 1.0)
            };
            let drop = load * share;
            world.elevation[idx] += drop;
            world.sediment[idx] += drop;
            budget.deposited += drop as f64;

            self.carry[idx] = 0.0;
            self.carry[r as usize] += load - drop;
        }

        // Anything still in transit ran off the edge of the world.
        for idx in 0..world.len() {
            budget.exported += self.carry[idx] as f64;
            self.carry[idx] = 0.0;
        }
    }

    /// Hillslope creep: nudge each cell toward the average of its neighbours.
    fn diffuse(&mut self, world: &mut World, dt: f32) {
        let p = &self.params;
        let sea = self.sea_level;
        if p.diffusion <= 0.0 {
            return;
        }
        let elevation = &world.elevation;

        self.scratch
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, out)| {
                let mut sum = 0.0;
                let mut n = 0.0;
                for (nidx, dx, dy) in world.neighbors8(idx) {
                    let weight = if dx != 0 && dy != 0 { 0.5 } else { 1.0 };
                    sum += elevation[nidx] * weight;
                    n += weight;
                }
                if n == 0.0 {
                    *out = elevation[idx];
                    return;
                }
                let mean = sum / n;
                let rate = if elevation[idx] > sea {
                    p.diffusion
                } else {
                    p.diffusion * p.submarine_diffusion_scale
                };
                *out = elevation[idx] + (mean - elevation[idx]) * (rate * dt).clamp(0.0, 1.0);
            });

        world.elevation.copy_from_slice(&self.scratch);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hydrology::{Hydrology, HydrologyParams};

    /// A ramp draining south, with a bit of texture so the flow concentrates.
    fn ramp(w: usize, h: usize) -> World {
        let mut world = World::new(w, h);
        for i in 0..world.len() {
            let (x, y) = world.coords(i);
            world.elevation[i] =
                6.0 - y as f32 * 0.12 + ((x as f32 * 0.7).sin() + (y as f32 * 0.4).cos()) * 0.15;
        }
        world
    }

    fn route(world: &mut World) -> FlowNetwork {
        let hydro = Hydrology::new(HydrologyParams::default(), 0.0, 5, world.width);
        let mut net = FlowNetwork::new(world.len());
        hydro.rainfall(world);
        hydro.route(world, &mut net);
        net
    }

    #[test]
    fn rivers_cut_down_into_a_slope() {
        let mut world = ramp(64, 48);
        let before = world.clone();
        let mut erosion = Erosion::new(ErosionParams::default(), 0.0, world.len());
        for _ in 0..20 {
            let net = route(&mut world);
            erosion.apply(&mut world, &net, 1.0);
        }
        let lowered = (0..world.len())
            .filter(|i| world.elevation[*i] < before.elevation[*i] - 1e-3)
            .count();
        assert!(
            lowered > world.len() / 10,
            "erosion barely did anything: {lowered} cells lowered"
        );
    }

    #[test]
    fn bigger_rivers_cut_deeper() {
        let mut world = ramp(64, 48);
        let mut erosion = Erosion::new(ErosionParams::default(), 0.0, world.len());
        // Let a drainage network develop first.
        for _ in 0..20 {
            let net = route(&mut world);
            erosion.apply(&mut world, &net, 1.0);
        }

        // Then measure the incision itself, not the net elevation change -
        // the trunk of a river is also where sediment gets dropped, so the two
        // are different questions.
        let net = route(&mut world);
        erosion.apply(&mut world, &net, 1.0);

        let mut samples: Vec<(f32, f32)> = (0..world.len())
            .filter(|i| net.receiver[*i] >= 0 && world.water[*i] == 0.0)
            .map(|i| (world.flow[i], erosion.cut[i]))
            .collect();
        samples.sort_by(|a, b| a.0.total_cmp(&b.0));

        let q = samples.len() / 4;
        let mean = |s: &[(f32, f32)]| s.iter().map(|v| v.1).sum::<f32>() / s.len() as f32;
        let small = mean(&samples[..q]);
        let large = mean(&samples[samples.len() - q..]);
        assert!(
            large > small * 2.0,
            "high discharge cut {large} should clearly beat low discharge {small}"
        );
    }

    #[test]
    fn sediment_is_conserved() {
        let mut world = ramp(48, 32);
        let mut erosion = Erosion::new(
            ErosionParams {
                diffusion: 0.0, // creep moves material too; isolate the rivers
                ..Default::default()
            },
            0.0,
            world.len(),
        );
        for _ in 0..10 {
            let net = route(&mut world);
            let budget = erosion.apply(&mut world, &net, 1.0);
            let accounted = budget.deposited + budget.exported;
            assert!(
                (budget.eroded - accounted).abs() <= budget.eroded * 1e-3 + 1e-6,
                "{} eroded but {} accounted for",
                budget.eroded,
                accounted
            );
        }
    }

    #[test]
    fn lakes_and_river_mouths_collect_sediment() {
        // A ramp running into a basin that has no outlet to the sea.
        let mut world = World::new(48, 48);
        for i in 0..world.len() {
            let (_, y) = world.coords(i);
            world.elevation[i] = 6.0 - y as f32 * 0.1;
        }
        for y in 30..38 {
            for x in 18..30 {
                let i = world.idx(x, y);
                world.elevation[i] = 1.0;
            }
        }
        let mut erosion = Erosion::new(ErosionParams::default(), 0.0, world.len());
        for _ in 0..30 {
            let net = route(&mut world);
            erosion.apply(&mut world, &net, 1.0);
        }
        let basin: f32 = (30..38)
            .flat_map(|y| (18..30).map(move |x| (x, y)))
            .map(|(x, y)| world.sediment[world.idx(x, y)])
            .sum();
        assert!(basin > 0.0, "the basin should be silting up");
    }

    #[test]
    fn the_sea_floor_is_left_alone_by_rivers() {
        let mut world = World::new(48, 32);
        world.elevation.fill(-3.0);
        for y in 0..8 {
            for x in 0..world.width {
                let i = world.idx(x, y);
                world.elevation[i] = 4.0 - y as f32 * 0.4;
            }
        }
        let before = world.clone();
        let mut erosion = Erosion::new(
            ErosionParams {
                diffusion: 0.0,
                ..Default::default()
            },
            0.0,
            world.len(),
        );
        for _ in 0..15 {
            let net = route(&mut world);
            erosion.apply(&mut world, &net, 1.0);
        }
        for y in 12..32 {
            for x in 0..world.width {
                let i = world.idx(x, y);
                assert!(
                    world.elevation[i] >= before.elevation[i] - 1e-6,
                    "deep sea floor was incised at ({x},{y})"
                );
            }
        }
    }

    #[test]
    fn diffusion_rounds_off_a_spike() {
        let mut world = World::new(32, 32);
        world.elevation.fill(1.0);
        let peak = world.idx(16, 16);
        world.elevation[peak] = 20.0;
        // An explicit rate, so this tests the diffusion itself rather than
        // whatever the default happens to be tuned to today.
        let mut erosion = Erosion::new(
            ErosionParams {
                diffusion: 0.2,
                ..Default::default()
            },
            0.0,
            world.len(),
        );
        for _ in 0..20 {
            erosion.diffuse(&mut world, 1.0);
        }
        assert!(world.elevation[peak] < 8.0, "spike survived: {}", world.elevation[peak]);
        assert!(world.elevation[world.idx(15, 16)] > 1.0, "material went nowhere");
    }

    #[test]
    fn stays_stable_under_long_runs() {
        let mut world = ramp(64, 48);
        let mut erosion = Erosion::new(ErosionParams::default(), 0.0, world.len());
        for _ in 0..200 {
            let net = route(&mut world);
            erosion.apply(&mut world, &net, 1.0);
        }
        for (i, e) in world.elevation.iter().enumerate() {
            assert!(e.is_finite(), "cell {i} went non-finite");
            assert!((-50.0..50.0).contains(e), "cell {i} blew up to {e}");
        }
    }

    #[test]
    fn erosion_is_deterministic() {
        let mut a = ramp(48, 32);
        let mut b = a.clone();
        let mut ea = Erosion::new(ErosionParams::default(), 0.0, a.len());
        let mut eb = Erosion::new(ErosionParams::default(), 0.0, b.len());
        for _ in 0..10 {
            let na = route(&mut a);
            ea.apply(&mut a, &na, 1.0);
            let nb = route(&mut b);
            eb.apply(&mut b, &nb, 1.0);
        }
        assert_eq!(a.elevation, b.elevation);
    }
}
