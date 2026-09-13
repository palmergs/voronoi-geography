//! The fixed terrain grid.
//!
//! The grid never moves; plates move across it (design doc section 2.1).
//! All large fields are flat `Vec`s indexed by `y * width + x`.

use crate::math::{Vec2, wrap_delta};
use crate::plate::PlateId;

pub const DEFAULT_WIDTH: usize = 1024;
pub const DEFAULT_HEIGHT: usize = 512;

/// Sentinel for "no plate owns this cell yet".
pub const NO_PLATE: PlateId = PlateId::MAX;

#[derive(Clone, Debug)]
pub struct World {
    pub width: usize,
    pub height: usize,

    /// Metres-ish. Sea level is 0.0; the unit is whatever the tuning says.
    pub elevation: Vec<f32>,

    /// Hydrology fields, filled in by the hydrology/erosion stages.
    pub rainfall: Vec<f32>,
    /// Standing water depth (lakes) above `elevation`.
    pub water: Vec<f32>,
    /// Flow accumulation (upstream rainfall volume).
    pub flow: Vec<f32>,
    /// Sediment deposited by the erosion stage, kept for debugging/visuals.
    pub sediment: Vec<f32>,

    /// Current Voronoi ownership.
    pub plate_id: Vec<PlateId>,
    /// Geological memory: time since this cell's crust was (re)created.
    pub crust_age: Vec<f32>,
}

impl World {
    pub fn new(width: usize, height: usize) -> Self {
        let n = width * height;
        World {
            width,
            height,
            elevation: vec![0.0; n],
            rainfall: vec![0.0; n],
            water: vec![0.0; n],
            flow: vec![0.0; n],
            sediment: vec![0.0; n],
            plate_id: vec![NO_PLATE; n],
            crust_age: vec![0.0; n],
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.width * self.height
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }

    #[inline]
    pub fn coords(&self, idx: usize) -> (usize, usize) {
        (idx % self.width, idx / self.width)
    }

    /// East/west wrap.
    #[inline]
    pub fn wrap_x(&self, x: i64) -> usize {
        x.rem_euclid(self.width as i64) as usize
    }

    /// North/south walls: out-of-range rows simply do not exist.
    #[inline]
    pub fn clamp_y(&self, y: i64) -> usize {
        y.clamp(0, self.height as i64 - 1) as usize
    }

    /// Index of a neighbour, wrapping in x and rejecting rows off the walls.
    #[inline]
    pub fn neighbor(&self, x: usize, y: usize, dx: i64, dy: i64) -> Option<usize> {
        let ny = y as i64 + dy;
        if ny < 0 || ny >= self.height as i64 {
            return None;
        }
        let nx = self.wrap_x(x as i64 + dx);
        Some(self.idx(nx, ny as usize))
    }

    /// Shortest vector from `(bx, by)` to `(ax, ay)` on the cylinder.
    #[inline]
    pub fn delta(&self, ax: f32, ay: f32, bx: f32, by: f32) -> Vec2 {
        Vec2::new(wrap_delta(ax, bx, self.width as f32), ay - by)
    }

    /// Squared distance between two points on the cylinder.
    #[inline]
    pub fn dist2(&self, ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
        self.delta(ax, ay, bx, by).length_squared()
    }

    /// Clamp a position to the world: wrap x, keep y inside the walls.
    pub fn contain(&self, p: Vec2) -> Vec2 {
        let w = self.width as f32;
        Vec2::new(p.x.rem_euclid(w), p.y.clamp(0.0, self.height as f32 - 1.0))
    }

    /// The 8 neighbours of a cell, wrapping in x and clipped at the walls.
    pub fn neighbors8(&self, idx: usize) -> Neighbors8<'_> {
        let (x, y) = self.coords(idx);
        Neighbors8 { world: self, x, y, i: 0 }
    }
}

/// Offsets and step costs for an 8-connected walk.
pub const OFFSETS8: [(i64, i64); 8] = [
    (1, 0),
    (1, 1),
    (0, 1),
    (-1, 1),
    (-1, 0),
    (-1, -1),
    (0, -1),
    (1, -1),
];

pub struct Neighbors8<'a> {
    world: &'a World,
    x: usize,
    y: usize,
    i: usize,
}

impl Iterator for Neighbors8<'_> {
    /// `(index, dx, dy)`
    type Item = (usize, i64, i64);

    fn next(&mut self) -> Option<Self::Item> {
        while self.i < OFFSETS8.len() {
            let (dx, dy) = OFFSETS8[self.i];
            self.i += 1;
            if let Some(n) = self.world.neighbor(self.x, self.y, dx, dy) {
                return Some((n, dx, dy));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn w() -> World {
        World::new(16, 8)
    }

    #[test]
    fn indexing_round_trips() {
        let w = w();
        for y in 0..w.height {
            for x in 0..w.width {
                assert_eq!(w.coords(w.idx(x, y)), (x, y));
            }
        }
    }

    #[test]
    fn x_wraps_and_y_walls() {
        let w = w();
        assert_eq!(w.wrap_x(-1), 15);
        assert_eq!(w.wrap_x(16), 0);
        assert_eq!(w.neighbor(15, 3, 1, 0), Some(w.idx(0, 3)));
        assert_eq!(w.neighbor(0, 0, 0, -1), None);
        assert_eq!(w.neighbor(0, 7, 0, 1), None);
    }

    #[test]
    fn neighbors8_clips_at_the_walls() {
        let w = w();
        assert_eq!(w.neighbors8(w.idx(5, 0)).count(), 5);
        assert_eq!(w.neighbors8(w.idx(5, 7)).count(), 5);
        assert_eq!(w.neighbors8(w.idx(0, 3)).count(), 8);
    }

    #[test]
    fn delta_crosses_the_seam() {
        let w = w();
        let d = w.delta(0.0, 1.0, 15.0, 1.0);
        assert_eq!(d.x, 1.0);
    }

    #[test]
    fn contain_wraps_and_clamps() {
        let w = w();
        let p = w.contain(Vec2::new(-1.0, -5.0));
        assert_eq!(p.x, 15.0);
        assert_eq!(p.y, 0.0);
        let p = w.contain(Vec2::new(17.0, 100.0));
        assert_eq!(p.x, 1.0);
        assert_eq!(p.y, 7.0);
    }
}
