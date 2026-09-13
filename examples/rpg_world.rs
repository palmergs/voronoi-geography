//! Generating a world once and querying it, the way a game would.
//!
//! Run with:
//!
//! ```sh
//! cargo run --release --example rpg_world
//! ```

use voronoi_geology::boundary::BoundaryKind;
use voronoi_geology::plate::{CrustType, PlateParams};
use voronoi_geology::query::Sample;
use voronoi_geology::simulation::{Simulation, SimulationParams};

fn main() {
    // 1. Generate. Everything derives from the seed, so storing the seed and
    //    the parameters is enough to reproduce this world exactly.
    let mut world = Simulation::new(SimulationParams {
        seed: 20_260_912,
        plates: PlateParams {
            count: 32,
            continental_fraction: 0.4,
            ..Default::default()
        },
        ..Default::default()
    });
    world.run(400);

    // 2. From here the simulation is just a read-only model. Keep it for the
    //    life of the game and query it.
    let threshold = world.river_threshold();

    let place = world.sample(512, 256);
    println!("{}", describe(&place, threshold));

    // 3. Coordinates that may run off a region use `sample_at`: x wraps around
    //    the world, y returns None past the northern and southern walls.
    let (px, py) = (2, 40);
    for dy in -1..=1 {
        for dx in -1..=1 {
            match world.sample_at(px + dx, py + dy) {
                Some(s) => print!("{}", if s.is_land() { '#' } else { '~' }),
                None => print!(" "), // off the world
            }
        }
        println!();
    }

    // 4. Bulk work - building a tile map, picking sites - should iterate the
    //    fields directly rather than sampling cell by cell.
    let land = world
        .world
        .elevation
        .iter()
        .filter(|e| **e > world.sea_level())
        .count();
    println!(
        "\n{:.1}% land, {} river cells",
        land as f32 / world.world.len() as f32 * 100.0,
        world.rivers(threshold).count()
    );

    // A settlement wants fresh water, workable land and no volcanoes.
    let site = world
        .rivers(threshold)
        .filter(|s| s.altitude() < 1.5 && s.rainfall > 0.4)
        .filter(|s| !matches!(s.boundary.map(|b| b.kind), Some(BoundaryKind::Convergent)))
        .max_by(|a, b| a.flow.total_cmp(&b.flow));

    match site {
        Some(s) => println!(
            "best river site: ({}, {}) at {:.0} m above sea level",
            s.x,
            s.y,
            s.altitude() * 1000.0
        ),
        None => println!("no suitable river site in this world"),
    }

    println!(
        "\nretains about {} MB of field data (peak while stepping is roughly \
         double that, for per-step scratch)",
        retained_mb(&world)
    );
}

/// Turn a sample into something a game could show a player.
fn describe(s: &Sample, river_threshold: f32) -> String {
    let ground = if !s.is_land() {
        match s.depth() {
            d if d > 4.0 => "deep ocean",
            d if d > 1.0 => "open sea",
            _ => "shallow water",
        }
    } else if s.is_lake() {
        "lake"
    } else {
        match s.altitude() {
            a if a > 4.5 => "high mountains",
            a if a > 2.0 => "mountains",
            a if a > 0.8 => "hills",
            _ => "lowland",
        }
    };

    let crust = match s.crust {
        CrustType::Continental => "continental crust",
        CrustType::Oceanic => "oceanic crust",
    };

    let setting = match s.boundary {
        Some(b) if b.distance < 12.0 => match b.kind {
            BoundaryKind::Convergent => ", in a zone of active mountain building",
            BoundaryKind::Divergent => ", where the crust is pulling apart",
            BoundaryKind::Transform => ", along a fault line",
        },
        _ => ", well inside a stable plate",
    };

    let water = if s.flow >= river_threshold {
        " A river runs through it."
    } else {
        ""
    };

    format!(
        "({}, {}): {ground} on {crust}{setting}.\n  \
         elevation {:+.2}, latitude {:.0} degrees, rainfall {:.2}, crust age {:.0}.{water}",
        s.x, s.y, s.elevation, s.latitude, s.rainfall, s.crust_age
    )
}

/// Field data a finished world holds on to, for budgeting.
///
/// This is what stays resident. Stepping allocates and frees a similar amount
/// again for scratch - the distance transform and the flood fill both need
/// working buffers - so peak usage during generation is roughly twice this.
/// Measured peak for the default 1024x512 world is around 70 MB.
fn retained_mb(sim: &Simulation) -> usize {
    let cells = sim.world.len();
    // 6 f32 fields + u16 plate ids, the flow network, the boundary distance
    // field, the boundary warp and the erosion scratch buffers.
    let per_cell = 6 * 4 + 2 + 12 + 8 + 8 + 12;
    (cells * per_cell) / (1024 * 1024)
}
