//! A game-oriented plate tectonics world simulation.
//!
//! See `docs/plate-tectonics-world-sim-project-design.md`. The short version:
//! a fixed 1024x512 cylindrical heightfield, with tectonic plates modelled as
//! moving Voronoi sites that drift across it. Boundary behaviour falls out of
//! the relative velocity of neighbouring plates, and every effect is applied
//! through a narrow distance falloff so features stay legible at this
//! resolution. Small, local, repeated - never large, global, instantaneous.

pub mod math;
pub mod noise;
pub mod plate;
pub mod voronoi;
pub mod world;
pub mod boundary;
pub mod tectonics;
pub mod hydrology;
pub mod erosion;
pub mod simulation;
pub mod render;
pub mod query;
