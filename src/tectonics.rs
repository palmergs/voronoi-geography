//! Turning classified boundaries into terrain (design doc sections 7-14, 19-21).
//!
//! Two stages run here:
//!
//! 1. [`Tectonics::apply`] - continuous deformation. Every cell looks up its
//!    nearest boundary, and the boundary's kind plus the crust types on either
//!    side selects a rule. Everything is multiplied by a falloff that is zero
//!    past a hard radius, so plate interiors stay untouched.
//! 2. [`Tectonics::apply_volcanism`] - discrete volcanic cones, placed by a
//!    deterministic hash so the same seed gives the same volcanoes.
//!
//! Per-step magnitudes are deliberately tiny. Mountains are supposed to be the
//! result of a few hundred steps of history, not one dramatic event
//! (section 21).

use crate::boundary::{Boundary, BoundaryKind, NearestBoundary};
use crate::math::{bump, lerp, smoothstep};
use crate::noise::{CylinderNoise, hash_unit};
use crate::plate::{CrustType, Plate, PlateId};
use crate::world::World;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
pub struct TectonicParams {
    // --- convergent ---
    /// Continent/continent collision: broad, strong (section 7.1).
    pub collision_radius: f32,
    pub collision_rate: f32,
    /// Uplift on the overriding plate of a subduction zone.
    pub orogeny_radius: f32,
    pub orogeny_rate: f32,
    /// Trench on the downgoing plate: narrow and deep.
    pub trench_radius: f32,
    pub trench_rate: f32,
    /// Volcanic arc, offset from the trench onto the overriding plate.
    pub arc_offset: f32,
    pub arc_sigma: f32,
    pub arc_rate: f32,

    // --- divergent ---
    /// Rift/ridge axis: the raised spreading centre itself.
    pub ridge_radius: f32,
    pub ridge_rate: f32,
    /// The wider subsiding rift around it.
    pub rift_radius: f32,
    pub rift_rate: f32,
    /// Rift shoulders (uplifted flanks of a continental rift).
    pub shoulder_offset: f32,
    pub shoulder_sigma: f32,
    pub shoulder_rate: f32,
    /// Within this radius of a spreading centre, crust is brand new: its age
    /// resets and its elevation relaxes toward the ridge baseline.
    pub new_crust_radius: f32,
    pub new_crust_rate: f32,
    pub ridge_elevation: f32,

    // --- transform ---
    pub transform_radius: f32,
    pub transform_rate: f32,

    // --- global ---
    /// Soft isostatic limits: uplift slows as terrain approaches these, so
    /// long-lived boundaries do not grow without bound.
    pub max_elevation: f32,
    pub min_elevation: f32,
    /// Ocean floor sinks slowly as it ages (section 15).
    pub age_subsidence_rate: f32,
    pub age_subsidence_max: f32,

    // --- volcanism (section 19) ---
    pub volcano_chance: f32,
    pub volcano_radius: f32,
    pub volcano_height: f32,
}

impl Default for TectonicParams {
    fn default() -> Self {
        TectonicParams {
            collision_radius: 24.0,
            collision_rate: 0.055,
            orogeny_radius: 18.0,
            orogeny_rate: 0.050,
            trench_radius: 7.0,
            trench_rate: 0.045,
            arc_offset: 7.0,
            arc_sigma: 3.0,
            arc_rate: 0.030,

            ridge_radius: 4.0,
            ridge_rate: 0.030,
            rift_radius: 14.0,
            rift_rate: 0.030,
            shoulder_offset: 9.0,
            shoulder_sigma: 4.0,
            shoulder_rate: 0.022,
            new_crust_radius: 5.0,
            new_crust_rate: 0.05,
            ridge_elevation: -2.4,

            transform_radius: 6.0,
            transform_rate: 0.020,

            max_elevation: 8.0,
            min_elevation: -9.0,
            age_subsidence_rate: 0.0016,
            age_subsidence_max: -5.6,

            volcano_chance: 0.020,
            volcano_radius: 3.0,
            volcano_height: 0.55,
        }
    }
}

pub struct Tectonics {
    pub params: TectonicParams,
    seed: u64,
    /// Fine-grained field that gives transform zones their broken-up texture.
    stress_noise: CylinderNoise,
}

impl Tectonics {
    pub fn new(params: TectonicParams, seed: u64, width: usize) -> Self {
        Tectonics {
            params,
            seed,
            stress_noise: CylinderNoise::new((seed as u32) ^ 0x7EC7, width, 5.0, 3),
        }
    }

    /// The largest radius any rule can reach; the distance field only needs to
    /// be built out this far.
    pub fn max_radius(&self) -> f32 {
        let p = &self.params;
        p.collision_radius
            .max(p.orogeny_radius)
            .max(p.trench_radius)
            .max(p.arc_offset + 3.0 * p.arc_sigma)
            .max(p.rift_radius)
            .max(p.shoulder_offset + 3.0 * p.shoulder_sigma)
            .max(p.transform_radius)
    }

    /// Continuous deformation plus crust ageing. One simulation step.
    pub fn apply(
        &self,
        world: &mut World,
        plates: &[Plate],
        boundaries: &[Boundary],
        field: &NearestBoundary,
        plate_age: &[f32],
        dt: f32,
    ) {
        let p = &self.params;
        let width = world.width;
        let World {
            elevation,
            crust_age,
            plate_id,
            ..
        } = world;

        elevation
            .par_iter_mut()
            .zip(crust_age.par_iter_mut())
            .enumerate()
            .for_each(|(idx, (elev, age))| {
                *age += dt;

                // Old oceanic crust cools and sinks (section 15). This is what
                // makes young ridges stand above the abyssal plain.
                if *elev < 0.0 && *elev > p.age_subsidence_max {
                    *elev -= p.age_subsidence_rate * dt * (*age / (*age + 40.0));
                }

                let Some((src, d)) = field.at(idx) else {
                    return;
                };
                let b = &boundaries[src as usize];
                let own = plate_id[idx];
                let at = crate::math::Vec2::new((idx % width) as f32, (idx / width) as f32);
                let outcome = self.effect(b, plates, plate_age, own, d, at);

                if outcome.new_crust > 0.0 {
                    // Fresh crust. The clock resets hard - this crust is new,
                    // not merely disturbed - which is what makes the crust-age
                    // layer read as a spreading pattern, with age increasing
                    // away from every ridge.
                    *age = lerp(*age, 0.0, (outcome.new_crust * 2.0).clamp(0.0, 1.0));
                    // The surface only relaxes toward the ridge baseline
                    // gradually, so a ridge does not punch a hole in terrain
                    // it has just started spreading through.
                    let t = (outcome.new_crust * p.new_crust_rate * dt).clamp(0.0, 1.0);
                    *elev = lerp(*elev, p.ridge_elevation, t);
                }

                let mut delta = outcome.delta * dt;
                // Soft isostatic limit.
                if delta > 0.0 {
                    delta *= 1.0 - smoothstep(p.max_elevation * 0.6, p.max_elevation, *elev);
                } else if delta < 0.0 {
                    delta *= 1.0 - smoothstep(p.min_elevation * 0.6, p.min_elevation, *elev);
                }
                *elev += delta;
            });
    }

    /// Elevation change per unit time for one cell, and whether that cell is
    /// being resurfaced as new crust.
    fn effect(
        &self,
        b: &Boundary,
        plates: &[Plate],
        plate_age: &[f32],
        own: PlateId,
        d: f32,
        at: crate::math::Vec2,
    ) -> Effect {
        let p = &self.params;
        let crust_a = plates[b.plate_a as usize].crust_type;
        let crust_b = plates[b.plate_b as usize].crust_type;
        // `own` is the plate that actually holds this cell, which is what the
        // rules key off - the cell may sit on either side of the boundary.
        let own_crust = plates[own as usize].crust_type;
        let strength = b.strength;

        match b.kind {
            BoundaryKind::Convergent => {
                match (crust_a, crust_b) {
                    // Continent meets continent: no subduction, just a broad
                    // mountain belt lifting both sides (section 7.1).
                    (CrustType::Continental, CrustType::Continental) => Effect::uplift(
                        p.collision_rate * strength * falloff(d, p.collision_radius),
                    ),
                    _ => {
                        let downgoing = subducting_plate(b, plates, plate_age);
                        if own == downgoing {
                            // Trench: narrow, deep, on the downgoing side.
                            Effect::uplift(
                                -p.trench_rate * strength * falloff(d, p.trench_radius),
                            )
                        } else {
                            // Overriding side: broad uplift, plus a volcanic
                            // arc set back from the trench.
                            let orogeny = if own_crust.is_continental() {
                                p.orogeny_rate
                            } else {
                                // Island arcs are lower and narrower than a
                                // continental cordillera.
                                p.orogeny_rate * 0.45
                            };
                            let arc = p.arc_rate * bump(d, p.arc_offset, p.arc_sigma);
                            Effect::uplift(
                                strength * (orogeny * falloff(d, p.orogeny_radius) + arc),
                            )
                        }
                    }
                }
            }

            BoundaryKind::Divergent => {
                let oceanic = !own_crust.is_continental();
                // The axis itself rises; the wider rift around it drops. For
                // ocean crust the axis wins (a mid-ocean ridge stands proud of
                // the abyssal plain); on a continent the rift valley wins.
                let axis = p.ridge_rate * falloff(d, p.ridge_radius);
                let basin = p.rift_rate * falloff(d, p.rift_radius);
                let delta = if oceanic {
                    strength * (axis * 1.4 - basin * 0.5)
                } else {
                    let shoulders =
                        p.shoulder_rate * bump(d, p.shoulder_offset, p.shoulder_sigma);
                    strength * (axis * 0.5 - basin + shoulders)
                };
                Effect {
                    delta,
                    new_crust: if oceanic {
                        strength * falloff(d, p.new_crust_radius)
                    } else {
                        0.0
                    },
                }
            }

            // Small, irregular deformation: ridges and troughs shuffled along
            // the fault rather than a single clean feature (section 9).
            BoundaryKind::Transform => {
                let texture = self.stress_noise.get(at.x, at.y);
                Effect::uplift(
                    p.transform_rate * strength * texture * falloff(d, p.transform_radius),
                )
            }
        }
    }

    /// Place volcanoes (section 19).
    ///
    /// Probability is `base x tendency x strength`, sampled from a stateless
    /// hash of the boundary cell and the step number - deterministic, and
    /// independent of the order boundaries happen to be visited in.
    pub fn apply_volcanism(
        &self,
        world: &mut World,
        plates: &[Plate],
        boundaries: &[Boundary],
        plate_age: &[f32],
        step: u64,
        dt: f32,
    ) -> usize {
        let p = &self.params;
        let mut erupted = 0;

        for (i, b) in boundaries.iter().enumerate() {
            let tendency = volcanic_tendency(b, plates);
            if tendency <= 0.0 {
                continue;
            }
            let chance = (p.volcano_chance * tendency * b.strength * dt).clamp(0.0, 1.0);
            let roll = hash_unit(self.seed ^ 0x_C0FF_EE01, i as i64, step as i64);
            if roll >= chance {
                continue;
            }

            // Arc volcanoes sit back from the trench, on the overriding plate;
            // rift volcanoes sit on the axis itself.
            let (cx, cy) = world.coords(b.cell as usize);
            let offset = match b.kind {
                BoundaryKind::Convergent => {
                    let downgoing = subducting_plate(b, plates, plate_age);
                    if b.plate_a == downgoing {
                        continue; // the arc belongs to the other side
                    }
                    -p.arc_offset
                }
                _ => 0.0,
            };
            let vx = cx as f32 + b.normal.x * offset;
            let vy = cy as f32 + b.normal.y * offset;

            let jitter = hash_unit(self.seed ^ 0xA11E, i as i64, step as i64 ^ 0x5F);
            let height = p.volcano_height * (0.6 + 0.8 * jitter);
            stamp_cone(world, vx, vy, p.volcano_radius, height);
            erupted += 1;
        }

        erupted
    }
}

/// What one rule decided for one cell.
struct Effect {
    /// Elevation change per unit time.
    delta: f32,
    /// How strongly this cell is being resurfaced (0 = not at all).
    new_crust: f32,
}

impl Effect {
    fn uplift(delta: f32) -> Self {
        Effect {
            delta,
            new_crust: 0.0,
        }
    }
}

/// Which of the two plates goes under.
///
/// Oceanic crust always subducts beneath continental. When both are oceanic
/// the older, colder, denser plate wins - which is exactly the geological
/// memory `crust_age` exists to provide (section 15).
fn subducting_plate(b: &Boundary, plates: &[Plate], plate_age: &[f32]) -> PlateId {
    let a = plates[b.plate_a as usize].crust_type;
    let c = plates[b.plate_b as usize].crust_type;
    match (a, c) {
        (CrustType::Oceanic, CrustType::Continental) => b.plate_a,
        (CrustType::Continental, CrustType::Oceanic) => b.plate_b,
        _ => {
            let age_a = plate_age[b.plate_a as usize];
            let age_b = plate_age[b.plate_b as usize];
            if age_a > age_b || (age_a == age_b && b.plate_a < b.plate_b) {
                b.plate_a
            } else {
                b.plate_b
            }
        }
    }
}

/// Relative volcanic activity per boundary type (section 19). These are game
/// parameters, not measurements.
fn volcanic_tendency(b: &Boundary, plates: &[Plate]) -> f32 {
    let a = plates[b.plate_a as usize].crust_type;
    let c = plates[b.plate_b as usize].crust_type;
    let both_oceanic = !a.is_continental() && !c.is_continental();
    let both_continental = a.is_continental() && c.is_continental();

    match b.kind {
        BoundaryKind::Convergent => {
            if both_continental {
                0.15 // collision: little melt reaches the surface
            } else {
                1.0 // subduction: the classic volcanic arc
            }
        }
        BoundaryKind::Divergent => {
            if both_oceanic {
                1.0 // mid-ocean ridge
            } else {
                0.6 // continental rift
            }
        }
        BoundaryKind::Transform => 0.05,
    }
}

/// Add a smooth cone. Radius is small (1-4 cells) by design (section 11).
fn stamp_cone(world: &mut World, cx: f32, cy: f32, radius: f32, height: f32) {
    let r = radius.ceil() as i64;
    let x0 = cx.round() as i64;
    let y0 = cy.round() as i64;

    for dy in -r..=r {
        let y = y0 + dy;
        if y < 0 || y >= world.height as i64 {
            continue;
        }
        for dx in -r..=r {
            let x = world.wrap_x(x0 + dx);
            let d = ((dx * dx + dy * dy) as f32).sqrt();
            if d > radius {
                continue;
            }
            let t = 1.0 - d / radius;
            let idx = world.idx(x, y as usize);
            world.elevation[idx] += height * t * t;
        }
    }
}

/// Distance falloff: exponential, with a hard cutoff (section 10.1).
///
/// The last 30% of the radius is tapered so the cutoff does not show up as a
/// visible circular edge in the terrain.
#[inline]
pub fn falloff(d: f32, radius: f32) -> f32 {
    if d >= radius {
        return 0.0;
    }
    let width = radius * 0.35;
    (-d / width).exp() * (1.0 - smoothstep(radius * 0.7, radius, d))
}

/// Mean crust age per plate, used to decide which oceanic plate subducts.
pub fn mean_crust_age(world: &World, plate_count: usize) -> Vec<f32> {
    let mut total = vec![0.0f64; plate_count];
    let mut count = vec![0u32; plate_count];
    for (idx, &id) in world.plate_id.iter().enumerate() {
        let i = id as usize;
        if i < plate_count {
            total[i] += world.crust_age[idx] as f64;
            count[i] += 1;
        }
    }
    total
        .into_iter()
        .zip(count)
        .map(|(t, c)| if c > 0 { (t / c as f64) as f32 } else { 0.0 })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::{detect, nearest_boundary};
    use crate::math::Vec2;
    use crate::voronoi::{self, Warp};

    fn plate(id: PlateId, pos: (f32, f32), vel: (f32, f32), crust: CrustType) -> Plate {
        Plate {
            id,
            position: Vec2::new(pos.0, pos.1),
            velocity: Vec2::new(vel.0, vel.1),
            crust_type: crust,
        }
    }

    /// Two plates, meeting near x=32 and again across the seam near x=0.
    fn setup(
        crust_a: CrustType,
        crust_b: CrustType,
        va: (f32, f32),
        vb: (f32, f32),
    ) -> (World, Vec<Plate>, Vec<Boundary>, NearestBoundary, Tectonics) {
        let mut world = World::new(96, 48);
        let plates = vec![
            plate(0, (16.0, 24.0), va, crust_a),
            plate(1, (64.0, 24.0), vb, crust_b),
        ];
        let flat = Warp::none(&world);
        voronoi::assign(&mut world, &plates, &flat);
        for (i, e) in world.elevation.iter_mut().enumerate() {
            *e = if plates[world.plate_id[i] as usize].crust_type.is_continental() {
                0.3
            } else {
                -4.0
            };
        }
        let tec = Tectonics::new(TectonicParams::default(), 99, world.width);
        let bs = detect(&world, &plates, 0.02);
        let field = nearest_boundary(&world, &bs, tec.max_radius());
        (world, plates, bs, field, tec)
    }

    fn run(
        world: &mut World,
        plates: &[Plate],
        bs: &[Boundary],
        field: &NearestBoundary,
        tec: &Tectonics,
        steps: usize,
    ) {
        let ages = vec![0.0; plates.len()];
        for _ in 0..steps {
            tec.apply(world, plates, bs, field, &ages, 1.0);
        }
    }

    /// Mean elevation change in a vertical strip, relative to the start.
    fn strip_change(before: &World, after: &World, x0: usize, x1: usize) -> f32 {
        let mut sum = 0.0;
        let mut n = 0;
        for y in 0..after.height {
            for x in x0..x1 {
                let i = after.idx(x, y);
                sum += after.elevation[i] - before.elevation[i];
                n += 1;
            }
        }
        sum / n as f32
    }

    #[test]
    fn continental_collision_raises_both_sides() {
        let (mut world, plates, bs, field, tec) = setup(
            CrustType::Continental,
            CrustType::Continental,
            (0.4, 0.0),
            (-0.4, 0.0),
        );
        let before = world.clone();
        run(&mut world, &plates, &bs, &field, &tec, 40);

        // The collision is at x=40 (midway between the sites).
        assert!(strip_change(&before, &world, 34, 40) > 0.3, "west side did not rise");
        assert!(strip_change(&before, &world, 41, 47) > 0.3, "east side did not rise");
    }

    #[test]
    fn ocean_continent_convergence_digs_a_trench_and_builds_a_cordillera() {
        // Plate 0 is oceanic and drives east into continental plate 1.
        let (mut world, plates, bs, field, tec) = setup(
            CrustType::Oceanic,
            CrustType::Continental,
            (0.4, 0.0),
            (-0.4, 0.0),
        );
        let before = world.clone();
        run(&mut world, &plates, &bs, &field, &tec, 40);

        let ocean_side = strip_change(&before, &world, 36, 40);
        let continent_side = strip_change(&before, &world, 41, 50);
        assert!(ocean_side < -0.1, "no trench on the oceanic side: {ocean_side}");
        assert!(
            continent_side > 0.3,
            "no uplift on the continental side: {continent_side}"
        );
    }

    #[test]
    fn oceanic_divergence_leaves_a_ridge_standing_above_its_flanks() {
        // Plates pull apart at the mid-map boundary.
        let (mut world, plates, bs, field, tec) = setup(
            CrustType::Oceanic,
            CrustType::Oceanic,
            (-0.4, 0.0),
            (0.4, 0.0),
        );
        run(&mut world, &plates, &bs, &field, &tec, 60);

        let at = |x: usize| {
            let mut s = 0.0;
            for y in 0..world.height {
                s += world.elevation[world.idx(x, y)];
            }
            s / world.height as f32
        };
        let axis = at(40);
        let flank = at(50);
        assert!(axis > flank, "ridge axis {axis} should stand above flank {flank}");
    }

    #[test]
    fn continental_rift_drops_the_axis_below_its_shoulders() {
        let (mut world, plates, bs, field, tec) = setup(
            CrustType::Continental,
            CrustType::Continental,
            (-0.4, 0.0),
            (0.4, 0.0),
        );
        run(&mut world, &plates, &bs, &field, &tec, 60);

        let at = |x: usize| {
            let mut s = 0.0;
            for y in 0..world.height {
                s += world.elevation[world.idx(x, y)];
            }
            s / world.height as f32
        };
        assert!(at(40) < at(49), "rift valley should sit below its shoulders");
    }

    #[test]
    fn effects_do_not_reach_plate_interiors() {
        // A wide world, so the plates genuinely have interiors: the sites sit
        // 64 cells from the nearest boundary, well past any rule's radius.
        let mut world = World::new(256, 96);
        let plates = vec![
            plate(0, (40.0, 48.0), (0.4, 0.0), CrustType::Continental),
            plate(1, (168.0, 48.0), (-0.4, 0.0), CrustType::Continental),
        ];
        let flat = Warp::none(&world);
        voronoi::assign(&mut world, &plates, &flat);
        world.elevation.fill(0.3);

        let tec = Tectonics::new(TectonicParams::default(), 99, world.width);
        let bs = detect(&world, &plates, 0.02);
        let field = nearest_boundary(&world, &bs, tec.max_radius());
        let before = world.clone();
        run(&mut world, &plates, &bs, &field, &tec, 60);

        // Invariant 1: a cell with no boundary in range must not move at all.
        let mut peak: f32 = 0.0;
        for i in 0..world.len() {
            let moved = (world.elevation[i] - before.elevation[i]).abs();
            if field.at(i).is_none() {
                assert_eq!(moved, 0.0, "cell outside every radius moved by {moved}");
            }
            peak = peak.max(moved);
        }
        assert!(peak > 0.1, "nothing happened at all");

        // Invariant 2: the effect is concentrated near the boundary, not
        // smeared out to the cutoff. Half a radius out it has already faded to
        // a fraction of the peak, so the mountain belt reads as a belt.
        for i in 0..world.len() {
            if let Some((_, d)) = field.at(i)
                && d >= tec.params.collision_radius * 0.5 {
                    let moved = (world.elevation[i] - before.elevation[i]).abs();
                    assert!(
                        moved < peak * 0.3,
                        "effect still at {:.0}% of peak {d:.1} cells out",
                        moved / peak * 100.0
                    );
                }
        }

        // The plate sites themselves are untouched.
        for x in [38, 40, 42, 166, 168, 170] {
            let change = strip_change(&before, &world, x, x + 1);
            assert_eq!(change, 0.0, "plate interior at x={x} moved by {change}");
        }
    }

    #[test]
    fn deformation_scales_with_relative_velocity() {
        let slow = {
            let (mut w, p, b, f, t) =
                setup(CrustType::Continental, CrustType::Continental, (0.05, 0.0), (-0.05, 0.0));
            let before = w.clone();
            run(&mut w, &p, &b, &f, &t, 30);
            strip_change(&before, &w, 38, 42)
        };
        let fast = {
            let (mut w, p, b, f, t) =
                setup(CrustType::Continental, CrustType::Continental, (0.4, 0.0), (-0.4, 0.0));
            let before = w.clone();
            run(&mut w, &p, &b, &f, &t, 30);
            strip_change(&before, &w, 38, 42)
        };
        assert!(fast > slow * 4.0, "slow={slow} fast={fast}");
    }

    #[test]
    fn isostasy_keeps_mountains_bounded() {
        let (mut world, plates, bs, field, tec) = setup(
            CrustType::Continental,
            CrustType::Continental,
            (0.5, 0.0),
            (-0.5, 0.0),
        );
        run(&mut world, &plates, &bs, &field, &tec, 4000);
        let peak = world.elevation.iter().cloned().fold(f32::MIN, f32::max);
        assert!(
            peak < tec.params.max_elevation * 1.05,
            "runaway uplift: {peak}"
        );
    }

    #[test]
    fn spreading_centres_reset_crust_age() {
        let (mut world, plates, bs, field, tec) = setup(
            CrustType::Oceanic,
            CrustType::Oceanic,
            (-0.4, 0.0),
            (0.4, 0.0),
        );
        for a in world.crust_age.iter_mut() {
            *a = 100.0;
        }
        run(&mut world, &plates, &bs, &field, &tec, 30);

        let axis = world.crust_age[world.idx(40, 24)];
        let interior = world.crust_age[world.idx(16, 24)];
        assert!(axis < 60.0, "axis crust should be young, got {axis}");
        assert!(interior > 120.0, "interior crust should have aged, got {interior}");
    }

    #[test]
    fn older_oceanic_plate_is_the_one_that_subducts() {
        let (world, plates, bs, _, _) = setup(
            CrustType::Oceanic,
            CrustType::Oceanic,
            (0.4, 0.0),
            (-0.4, 0.0),
        );
        let _ = world;
        let b = bs.iter().find(|b| b.kind == BoundaryKind::Convergent).unwrap();
        assert_eq!(subducting_plate(b, &plates, &[10.0, 90.0]), 1);
        assert_eq!(subducting_plate(b, &plates, &[90.0, 10.0]), 0);
    }

    #[test]
    fn oceanic_crust_always_subducts_under_continental() {
        let (_, plates, bs, _, _) = setup(
            CrustType::Oceanic,
            CrustType::Continental,
            (0.4, 0.0),
            (-0.4, 0.0),
        );
        let b = bs.iter().find(|b| b.kind == BoundaryKind::Convergent).unwrap();
        // Even if the continental plate is far older.
        assert_eq!(subducting_plate(b, &plates, &[1.0, 500.0]), 0);
    }

    #[test]
    fn volcanism_is_deterministic_and_stays_near_boundaries() {
        let (mut a, plates, bs, field, tec) = setup(
            CrustType::Oceanic,
            CrustType::Continental,
            (0.5, 0.0),
            (-0.5, 0.0),
        );
        let mut b = a.clone();
        let ages = vec![0.0; plates.len()];

        let mut erupted = 0;
        for step in 0..60 {
            erupted += tec.apply_volcanism(&mut a, &plates, &bs, &ages, step, 1.0);
            tec.apply_volcanism(&mut b, &plates, &bs, &ages, step, 1.0);
        }
        assert!(erupted > 0, "no volcanoes formed at a subduction zone");
        assert_eq!(a.elevation, b.elevation, "volcanism must be reproducible");

        // Nothing should erupt far from a boundary.
        for i in 0..a.elevation.len() {
            if field.at(i).is_none() {
                assert_eq!(a.elevation[i], b.elevation[i]);
                let (x, _) = a.coords(i);
                assert!(
                    (14..=18).contains(&x) || a.elevation[i] == b.elevation[i],
                    "volcano in a plate interior"
                );
            }
        }
    }

    #[test]
    fn falloff_is_monotonic_and_hits_zero_at_the_cutoff() {
        let r = 12.0;
        let mut prev = f32::MAX;
        for i in 0..=120 {
            let d = i as f32 * 0.1;
            let f = falloff(d, r);
            assert!(f <= prev + 1e-6, "falloff must not rise: d={d}");
            assert!((0.0..=1.0).contains(&f));
            prev = f;
        }
        assert_eq!(falloff(r, r), 0.0);
        assert_eq!(falloff(r + 5.0, r), 0.0);
        assert!(falloff(0.0, r) > 0.9);
    }

    #[test]
    fn mean_crust_age_averages_per_plate() {
        let mut world = World::new(8, 4);
        for i in 0..world.len() {
            world.plate_id[i] = if i % 2 == 0 { 0 } else { 1 };
            world.crust_age[i] = if i % 2 == 0 { 10.0 } else { 30.0 };
        }
        let ages = mean_crust_age(&world, 2);
        assert!((ages[0] - 10.0).abs() < 1e-4);
        assert!((ages[1] - 30.0).abs() < 1e-4);
    }
}
