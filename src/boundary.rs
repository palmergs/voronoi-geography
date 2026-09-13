//! Boundary detection and classification (design doc sections 6, 7, 8, 9).
//!
//! A boundary record is produced per (cell, neighbouring plate) pair along the
//! edge of a Voronoi region. Because a Voronoi edge is the perpendicular
//! bisector of its two sites, the boundary normal is simply the direction from
//! one site to the other - no gradient estimation needed, and it stays stable
//! as the diagram is rebuilt each step.
//!
//! Classification comes entirely from relative plate velocity, so changing a
//! plate's velocity changes the behaviour of all of its boundaries.

use crate::math::{Vec2, perp};
use crate::plate::{Plate, PlateId};
use crate::world::World;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryKind {
    Convergent,
    Divergent,
    Transform,
}

#[derive(Clone, Copy, Debug)]
pub struct Boundary {
    /// The cell this record sits on; it belongs to `plate_a`.
    pub cell: u32,
    pub plate_a: PlateId,
    /// The plate on the far side.
    pub plate_b: PlateId,
    /// Unit normal, pointing from a's site toward b's site.
    pub normal: Vec2,
    /// `dot(v_a - v_b, normal)`. Positive closes the gap, negative opens it.
    pub convergence: f32,
    /// Relative motion along the boundary itself.
    pub shear: f32,
    pub kind: BoundaryKind,
    /// Magnitude driving this boundary's effects (design doc section 12).
    pub strength: f32,
}

/// Marker for "no boundary within range" in a [`NearestBoundary`].
pub const NO_SOURCE: u32 = u32::MAX;

/// For every cell: the closest boundary record, and how far away it is.
///
/// Nearest-boundary attribution (rather than summing every boundary in range)
/// is what keeps features attributable: a player looking at a mountain range
/// can trace it to one boundary, which is the guiding principle in section 28.
pub struct NearestBoundary {
    pub dist: Vec<f32>,
    pub source: Vec<u32>,
    pub max_radius: f32,
}

impl NearestBoundary {
    /// The boundary affecting `idx`, if one is within range.
    #[inline]
    pub fn at(&self, idx: usize) -> Option<(u32, f32)> {
        let s = self.source[idx];
        if s == NO_SOURCE {
            None
        } else {
            Some((s, self.dist[idx]))
        }
    }
}

/// Find every boundary cell and classify it.
///
/// `transform_threshold` is the normal-motion magnitude below which a boundary
/// is treated as pure strike-slip.
pub fn detect(world: &World, plates: &[Plate], transform_threshold: f32) -> Vec<Boundary> {
    let mut out = Vec::new();
    const NEIGHBORS4: [(i64, i64); 4] = [(1, 0), (-1, 0), (0, 1), (0, -1)];

    for y in 0..world.height {
        for x in 0..world.width {
            let idx = world.idx(x, y);
            let a = world.plate_id[idx];
            let mut seen: [PlateId; 4] = [a; 4];
            let mut seen_len = 1;

            for (dx, dy) in NEIGHBORS4 {
                let Some(n) = world.neighbor(x, y, dx, dy) else {
                    continue;
                };
                let b = world.plate_id[n];
                if seen[..seen_len].contains(&b) {
                    continue;
                }
                seen[seen_len] = b;
                seen_len += 1;
                out.push(classify(world, plates, idx as u32, a, b, transform_threshold));
            }
        }
    }

    out
}

fn classify(
    world: &World,
    plates: &[Plate],
    cell: u32,
    a: PlateId,
    b: PlateId,
    transform_threshold: f32,
) -> Boundary {
    let pa = &plates[a as usize];
    let pb = &plates[b as usize];

    // On a cylinder two plates meet along *two* boundaries, and the correct
    // normal is opposite at each. Measuring both sites relative to this cell
    // picks the site images that are actually local to it, so each boundary
    // gets the normal its own geometry implies.
    let (cx, cy) = world.coords(cell as usize);
    let (cx, cy) = (cx as f32, cy as f32);
    let to_a = world.delta(pa.position.x, pa.position.y, cx, cy);
    let to_b = world.delta(pb.position.x, pb.position.y, cx, cy);
    let normal = (to_b - to_a).normalize_or_zero();
    let tangent = perp(normal);
    let relative = pa.velocity - pb.velocity;

    let convergence = relative.dot(normal);
    let shear = relative.dot(tangent);

    let (kind, strength) = if convergence > transform_threshold {
        (BoundaryKind::Convergent, convergence)
    } else if convergence < -transform_threshold {
        (BoundaryKind::Divergent, -convergence)
    } else {
        (BoundaryKind::Transform, shear.abs())
    };

    Boundary {
        cell,
        plate_a: a,
        plate_b: b,
        normal,
        convergence,
        shear,
        kind,
        strength,
    }
}

/// Distance to the nearest boundary cell, plus which boundary that was.
///
/// Exact Euclidean distance transform (Felzenszwalb & Huttenlocher): a
/// vertical sweep per column, then a lower-envelope pass per row. Both are
/// O(cells), and unlike a chamfer approximation the result is exact - which
/// matters because this distance is what every tectonic falloff is keyed to.
///
/// The x wrap is handled by running the row pass over three copies of the row
/// and keeping the middle one, so the nearest source may be reached the short
/// way around the cylinder.
pub fn nearest_boundary(world: &World, boundaries: &[Boundary], max_radius: f32) -> NearestBoundary {
    const INF: f32 = 1e20;
    let (w, h, n) = (world.width, world.height, world.len());

    // Seed cells. Several boundaries can share a cell at a triple junction;
    // the strongest one wins, so the visible feature is the dominant process.
    let mut seed = vec![NO_SOURCE; n];
    for (i, b) in boundaries.iter().enumerate() {
        let idx = b.cell as usize;
        let current = seed[idx];
        if current == NO_SOURCE || boundaries[current as usize].strength < b.strength {
            seed[idx] = i as u32;
        }
    }

    // Pass 1: nearest source within each column (walls, so no wrapping).
    let mut col_d2 = vec![INF; n];
    let mut col_src = vec![NO_SOURCE; n];
    for x in 0..w {
        let mut best_y = 0i64;
        let mut best_src = NO_SOURCE;
        for y in 0..h {
            let idx = y * w + x;
            if seed[idx] != NO_SOURCE {
                best_y = y as i64;
                best_src = seed[idx];
            }
            if best_src != NO_SOURCE {
                let d = (y as i64 - best_y) as f32;
                col_d2[idx] = d * d;
                col_src[idx] = best_src;
            }
        }
        best_src = NO_SOURCE;
        for y in (0..h).rev() {
            let idx = y * w + x;
            if seed[idx] != NO_SOURCE {
                best_y = y as i64;
                best_src = seed[idx];
            }
            if best_src != NO_SOURCE {
                let d = (best_y - y as i64) as f32;
                let d2 = d * d;
                if d2 < col_d2[idx] {
                    col_d2[idx] = d2;
                    col_src[idx] = best_src;
                }
            }
        }
    }

    // Pass 2: lower envelope along each row, over three tiled copies.
    let m = 3 * w;
    let mut dist = vec![f32::MAX; n];
    let mut source = vec![NO_SOURCE; n];

    dist.par_chunks_mut(w)
        .zip(source.par_chunks_mut(w))
        .enumerate()
        .for_each_init(
            || Scratch::new(m),
            |scratch, (y, (dist_row, source_row))| {
                for x in 0..m {
                    scratch.f[x] = col_d2[y * w + (x % w)];
                }
                lower_envelope(scratch);

                for x in 0..w {
                    let q = x + w; // the middle copy
                    let src = col_src[y * w + (scratch.arg[q] % w)];
                    if src == NO_SOURCE {
                        continue;
                    }
                    let d = scratch.d[q].max(0.0).sqrt();
                    // Hard cutoff: outside the radius nothing happens at all
                    // (design doc section 10.1).
                    if d <= max_radius {
                        dist_row[x] = d;
                        source_row[x] = src;
                    }
                }
            },
        );

    NearestBoundary {
        dist,
        source,
        max_radius,
    }
}

/// Reusable buffers for [`lower_envelope`], one per rayon worker.
struct Scratch {
    /// Input: squared distance to the nearest source in each column.
    f: Vec<f32>,
    /// Output: squared 2D distance.
    d: Vec<f32>,
    /// Output: the column that won, so the source index can be looked up.
    arg: Vec<usize>,
    /// Parabola positions.
    v: Vec<usize>,
    /// Parabola intersections.
    z: Vec<f32>,
}

impl Scratch {
    fn new(m: usize) -> Self {
        Scratch {
            f: vec![0.0; m],
            d: vec![0.0; m],
            arg: vec![0; m],
            v: vec![0; m],
            z: vec![0.0; m + 1],
        }
    }
}

/// 1D squared distance transform of a sampled function.
///
/// Computes `d[q] = min over p of ((q - p)^2 + f[p])`, plus the winning `p`.
fn lower_envelope(s: &mut Scratch) {
    let n = s.f.len();
    if n == 0 {
        return;
    }
    let f = &s.f;
    let (v, z) = (&mut s.v, &mut s.z);

    let mut k = 0usize;
    v[0] = 0;
    z[0] = f32::NEG_INFINITY;
    z[1] = f32::INFINITY;

    let intersect = |q: usize, p: usize| -> f32 {
        let (qf, pf) = (q as f32, p as f32);
        ((f[q] + qf * qf) - (f[p] + pf * pf)) / (2.0 * qf - 2.0 * pf)
    };

    for q in 1..n {
        let mut s_ = intersect(q, v[k]);
        while s_ <= z[k] && k > 0 {
            k -= 1;
            s_ = intersect(q, v[k]);
        }
        k += 1;
        v[k] = q;
        z[k] = s_;
        z[k + 1] = f32::INFINITY;
    }

    k = 0;
    for q in 0..n {
        while z[k + 1] < q as f32 {
            k += 1;
        }
        let p = v[k];
        let d = q as f32 - p as f32;
        s.d[q] = d * d + f[p];
        s.arg[q] = p;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plate::{CrustField, CrustType, PlateParams, generate};
    use crate::voronoi::{self, Warp};

    fn plate(id: PlateId, pos: (f32, f32), vel: (f32, f32)) -> Plate {
        Plate {
            id,
            position: Vec2::new(pos.0, pos.1),
            velocity: Vec2::new(vel.0, vel.1),
            crust_type: CrustType::Oceanic,
        }
    }

    /// Two sites side by side: the boundary runs vertically between them and
    /// the normal points east, so x velocities decide the classification.
    fn two_plate_world(va: (f32, f32), vb: (f32, f32)) -> (World, Vec<Plate>, Vec<Boundary>) {
        let mut world = World::new(64, 32);
        let plates = vec![plate(0, (16.0, 16.0), va), plate(1, (48.0, 16.0), vb)];
        let flat = Warp::none(&world);
        voronoi::assign(&mut world, &plates, &flat);
        let b = detect(&world, &plates, 0.05);
        (world, plates, b)
    }

    /// Kinds along the mid-map boundary (near x=32) and along the wrapped one
    /// (near x=0/64). On a cylinder these two are always opposites.
    fn kinds(world: &World, bs: &[Boundary]) -> (Vec<BoundaryKind>, Vec<BoundaryKind>) {
        let mut middle = Vec::new();
        let mut seam = Vec::new();
        for b in bs {
            let (x, _) = world.coords(b.cell as usize);
            if (24..40).contains(&x) {
                middle.push(b.kind);
            } else if !(8..=56).contains(&x) {
                seam.push(b.kind);
            }
        }
        assert!(!middle.is_empty() && !seam.is_empty());
        (middle, seam)
    }

    #[test]
    fn convergent_divergent_transform_come_from_relative_velocity() {
        // Plate 0 drives east, plate 1 drives west: they close in the middle
        // of the map and open up on the far side of the cylinder.
        let (world, _, bs) = two_plate_world((0.5, 0.0), (-0.5, 0.0));
        let (middle, seam) = kinds(&world, &bs);
        assert!(middle.iter().all(|k| *k == BoundaryKind::Convergent));
        assert!(seam.iter().all(|k| *k == BoundaryKind::Divergent));

        // Reverse both velocities and the two boundaries swap roles.
        let (world, _, bs) = two_plate_world((-0.5, 0.0), (0.5, 0.0));
        let (middle, seam) = kinds(&world, &bs);
        assert!(middle.iter().all(|k| *k == BoundaryKind::Divergent));
        assert!(seam.iter().all(|k| *k == BoundaryKind::Convergent));

        // Pure north/south shear against a north/south boundary: no normal
        // motion anywhere, so both boundaries are transform.
        let (_, _, sliding) = two_plate_world((0.0, 0.5), (0.0, -0.5));
        assert!(sliding.iter().all(|b| b.kind == BoundaryKind::Transform));
    }

    #[test]
    fn strength_scales_with_relative_motion() {
        let (_, _, slow) = two_plate_world((0.05, 0.0), (-0.05, 0.0));
        let (_, _, fast) = two_plate_world((0.5, 0.0), (-0.5, 0.0));
        let peak = |bs: &Vec<Boundary>| bs.iter().map(|b| b.strength).fold(0.0, f32::max);
        assert!(peak(&fast) > peak(&slow) * 5.0);
    }

    #[test]
    fn normals_point_from_a_toward_b() {
        let (world, plates, bs) = two_plate_world((0.5, 0.0), (-0.5, 0.0));
        for b in &bs {
            let (x, y) = world.coords(b.cell as usize);
            let pb = plates[b.plate_b as usize];
            let toward_b = world.delta(pb.position.x, pb.position.y, x as f32, y as f32);
            assert!(
                b.normal.dot(toward_b) > 0.0,
                "normal {:?} at ({x},{y}) does not face plate b at {:?}",
                b.normal,
                pb.position
            );
            assert!((b.normal.length() - 1.0).abs() < 1e-5);
        }
    }

    #[test]
    fn boundaries_are_found_on_both_sides_and_nowhere_else() {
        let (world, _, bs) = two_plate_world((0.5, 0.0), (-0.5, 0.0));
        assert!(!bs.is_empty());
        for b in &bs {
            let idx = b.cell as usize;
            assert_eq!(world.plate_id[idx], b.plate_a);
            let touches_other = world
                .neighbors8(idx)
                .any(|(n, _, _)| world.plate_id[n] == b.plate_b);
            assert!(touches_other, "boundary cell is not actually adjacent to plate b");
        }
    }

    #[test]
    fn boundaries_are_detected_across_the_seam() {
        let mut world = World::new(64, 16);
        // Sites at x=8 and x=40: one boundary near x=24, the other wraps at x=56.
        let plates = vec![plate(0, (8.0, 8.0), (0.2, 0.0)), plate(1, (40.0, 8.0), (-0.2, 0.0))];
        let flat = Warp::none(&world);
        voronoi::assign(&mut world, &plates, &flat);
        let bs = detect(&world, &plates, 0.05);
        let seam_side = bs
            .iter()
            .filter(|b| {
                let (x, _) = world.coords(b.cell as usize);
                x > 50
            })
            .count();
        assert!(seam_side > 0, "the wrapped boundary should be detected too");
    }

    #[test]
    fn distance_field_matches_brute_force() {
        let mut world = World::new(96, 48);
        let plates = generate(&world, &PlateParams { count: 6, ..Default::default() }, 21, &CrustField::new(&world, 21, 0.4));
        let flat = Warp::none(&world);
        voronoi::assign(&mut world, &plates, &flat);
        let bs = detect(&world, &plates, 0.05);
        let radius = 12.0;
        let field = nearest_boundary(&world, &bs, radius);

        for y in 0..world.height {
            for x in 0..world.width {
                let idx = world.idx(x, y);
                let mut truth = f32::MAX;
                for b in &bs {
                    let (bx, by) = world.coords(b.cell as usize);
                    let d = world
                        .delta(bx as f32, by as f32, x as f32, y as f32)
                        .length();
                    truth = truth.min(d);
                }
                match field.at(idx) {
                    // The transform is exact, so this is not a tolerance test.
                    Some((_, d)) => assert!(
                        (d - truth).abs() < 1e-3,
                        "distance wrong at ({x},{y}): {d} vs {truth}"
                    ),
                    None => assert!(
                        truth > radius,
                        "cell ({x},{y}) is {truth} away but was cut off at {radius}"
                    ),
                }
            }
        }
    }

    #[test]
    fn distance_field_respects_the_hard_cutoff() {
        let mut world = World::new(128, 64);
        let plates = generate(&world, &PlateParams { count: 4, ..Default::default() }, 3, &CrustField::new(&world, 3, 0.4));
        let flat = Warp::none(&world);
        voronoi::assign(&mut world, &plates, &flat);
        let bs = detect(&world, &plates, 0.05);
        let field = nearest_boundary(&world, &bs, 6.0);
        assert!(field.dist.iter().filter(|d| **d <= 6.0).count() > 0);
        for i in 0..world.len() {
            if let Some((s, d)) = field.at(i) {
                assert!(d <= 6.0);
                assert!((s as usize) < bs.len());
            }
        }
        // With only 4 plates and a 6 cell radius, most of the map is untouched.
        let affected = (0..world.len()).filter(|i| field.at(*i).is_some()).count();
        assert!(
            affected < world.len() / 2,
            "effects should stay local: {affected} of {} cells",
            world.len()
        );
    }
}
