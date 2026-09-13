# Plate Tectonics World Simulation

## 1. Project Overview

This project is a **game-oriented plate tectonics simulation** intended to generate believable procedural worlds.

The goal is **not** to reproduce real geophysics with high fidelity. Instead, the simulation should produce landscapes whose geological structure looks and feels as though it developed through a coherent tectonic history.

The simulation should favor:

- visually understandable geological rules
- deterministic procedural generation
- relatively simple algorithms
- good performance in Rust
- effects that remain visible at the target resolution
- geological history that can plausibly produce interesting terrain
- tunable parameters suitable for a game

The world is represented as a **cylindrical map**:

- East/West wraps continuously.
- North/South are literal world walls.
- The terrain is represented by a fixed `1024 × 512` heightfield.

---

# 2. Core World Representation

## 2.1 Terrain Grid

The terrain consists of:

```text
1024 × 512 = 524,288 points
```

Each point stores elevation and eventually additional geological/hydrological information.

The grid does **not** physically move when plates move.

Instead, the tectonic plates move across the fixed terrain grid. This keeps the simulation substantially easier to reason about and implement.

Basic indexing:

```rust
index = y * width + x
```

with:

```rust
width  = 1024
height = 512
```

East/West coordinates wrap:

```rust
x = x.rem_euclid(width)
```

North/South coordinates stop at the world boundaries.

---

# 3. Tectonic Plates

The world is divided into tectonic plates using a **Voronoi diagram**.

Each plate is represented by a moving Voronoi site.

A minimal plate might initially contain:

```rust
struct Plate {
    id: PlateId,
    position: Vec2,
    velocity: Vec2,
    crust_type: CrustType,
}
```

where:

```rust
enum CrustType {
    Oceanic,
    Continental,
}
```

The exact number of plates is still a tuning parameter, but approximately **20–50 plates** is a useful initial range.

---

# 4. Plate Movement

## 4.1 Linear Velocity

The initial model uses a simple velocity vector:

```rust
velocity = Vec2 {
    x: vx,
    y: vy,
}
```

This means the entire plate translates in the same direction at the same speed.

For example:

```text
velocity = (0.8, 0.2)
```

means the plate moves mostly eastward with a smaller northward component.

This is intentionally simpler than modeling spherical plate rotation.

---

## 4.2 Angular Velocity

Angular velocity was considered as a possible alternative.

With angular motion, a plate rotates around some point/axis. Consequently, different parts of the plate have different velocities:

```text
          ↑
       ↗  |  ↖
          |
←---------●---------→
          |
       ↘  |  ↙
          ↓
```

That is closer to how real tectonic plates behave on a sphere.

However, because this world is cylindrical and the goal is a game simulation rather than geological research, angular velocity is unnecessary for the first implementation.

### Decision

**Use linear `vx/vy` velocity for version 1.**

Angular/Euler-pole-style motion can be introduced later if the simpler model proves insufficient.

---

# 5. Plate Movement and Voronoi Recalculation

Each simulation iteration begins by moving the plate sites.

Conceptually:

```rust
plate.position += plate.velocity * dt;
```

The Voronoi diagram is then recalculated.

This causes the plate boundaries to move across the terrain.

The terrain grid itself remains fixed.

The resulting process is:

```text
Move plates
     ↓
Recalculate Voronoi ownership
     ↓
Identify plate boundaries
     ↓
Calculate relative plate motion
     ↓
Generate tectonic effects
     ↓
Apply effects to terrain
     ↓
Run volcanism
     ↓
Run erosion/hydrology
```

---

# 6. Boundary Classification

Boundary behavior should be derived from the **relative velocity of the two plates**, rather than manually assigning every boundary as convergent/divergent/transform.

For two adjacent plates:

```rust
relative_velocity = velocity_a - velocity_b;
```

The boundary has a normal direction.

The component of relative velocity along that normal determines the boundary behavior.

Conceptually:

```text
                 boundary
              ───────────────
                    ↑
                    │ normal
                    │
                    ↓
```

Calculate:

```text
normal_motion = dot(relative_velocity, boundary_normal)
```

Then:

```text
normal_motion > threshold
    → convergent

normal_motion < -threshold
    → divergent

otherwise
    → transform
```

This means changing a plate's velocity naturally changes the behavior of its boundaries.

---

# 7. Boundary Types

The simulation recognizes three primary boundary types.

## 7.1 Convergent

Two plates move toward one another.

The exact geological result depends heavily on crust type.

### Oceanic + Oceanic

Expected effects:

- deep-ish trench
- uplift near the opposing boundary
- volcanic island arc
- localized deformation

### Oceanic + Continental

Expected effects:

- strong trench on the oceanic side
- strong uplift on the continental side
- volcanic mountain chain
- broad deformation zone

### Continental + Continental

Expected effects:

- broad mountain-building region
- substantial uplift
- relatively wide deformation zone
- potentially large mountain chains

---

# 8. Divergent Boundaries

Two plates move apart.

Expected effects:

- depression/rift near the boundary
- localized uplift along the actual spreading center
- volcanism
- creation of young crust

Oceanic divergence should tend toward a **mid-ocean-ridge-like feature**.

Continental divergence should tend toward a **rift valley** with associated volcanism.

The boundary should remain visually distinct rather than producing a huge low-elevation area extending across the plates.

---

# 9. Transform Boundaries

Two plates move primarily parallel to their shared boundary.

Expected effects:

- relatively small elevation changes
- localized stress/fissuring
- small ridges/depressions
- potentially irregular terrain

Transform boundaries should generally be less dramatic than convergent boundaries.

---

# 10. Localized Boundary Effects

A major design requirement is that tectonic effects **must not spread too far into the interiors of plates**.

The terrain resolution is only:

```text
1024 × 512
```

so an overly broad falloff could cause an entire plate to become a distorted mountain/rift region.

The rules need to remain visually legible.

## 10.1 Falloff

Tectonic effects should use a distance-based falloff.

Conceptually:

```text
effect = strength × falloff(distance)
```

A useful starting model is exponential:

```text
falloff(d) = exp(-d / width)
```

combined with a hard cutoff.

For example:

```text
distance < radius
    apply effect

distance >= radius
    apply nothing
```

This gives us both:

- smooth boundaries
- a firm limit on how far the effect can propagate

---

# 11. Initial Falloff Scale

These are **game-tuning values**, not geological measurements.

A useful starting point is:

| Effect | Initial scale |
|---|---:|
| General boundary deformation | 5–20 cells |
| Mountain-building | ~12–30 cells |
| Volcanic features | ~1–4 cells |
| Transform stress | ~3–10 cells |

These should be exposed as parameters rather than hardcoded.

The most important principle is:

> A player should be able to look at a mountain/rift/trench and understand that it belongs to a particular plate boundary.

The model should avoid making entire plate interiors look as though they were directly affected by the boundary.

---

# 12. Relative Motion as Tectonic Strength

Boundary effects should ideally scale with the amount of relative plate motion.

For example:

```text
slow convergence
    → weak deformation

fast convergence
    → strong deformation
```

Rather than:

```text
every convergent boundary
    → exactly +100 elevation
```

This allows plate velocity to actually matter.

A generalized strength calculation could eventually look like:

```rust
let strength = relative_velocity.length() * boundary_factor;
```

with additional modifiers for crust type and boundary type.

---

# 13. Tectonic Effects

Rather than immediately modifying terrain for every geological process, tectonics should produce intermediate **effects**.

For example:

```rust
struct TectonicEffect {
    location: Vec2,
    radius: f32,
    uplift: f32,
    depression: f32,
    volcanism: f32,
    stress: f32,
}
```

The pipeline then becomes:

```text
Plate movement
      ↓
Boundary analysis
      ↓
TectonicEffect generation
      ↓
Terrain modification
```

This separation should make the system much easier to tune.

It also makes it possible to visualize/debug tectonic forces independently of the final terrain.

---

# 14. Crust Type

Every plate has a crust type:

```rust
enum CrustType {
    Oceanic,
    Continental,
}
```

Crust type affects the behavior of boundaries.

For example:

```text
Oceanic → Continental
    strong subduction behavior

Continental → Continental
    mountain building

Oceanic → Oceanic
    trenches + volcanic island arcs
```

The initial system does not need a complex physical density model.

Crust type is primarily a **rule selector**.

---

# 15. Crust Age

A potentially valuable addition is **crust age**.

Each terrain point can store:

```rust
crust_age: f32
```

or age could initially be associated with plates and later moved to individual terrain cells.

The important geological idea is:

```text
divergent boundary
       ↓
new crust
       ↓
crust moves away
       ↓
crust gets older
       ↓
eventually reaches a convergent boundary
       ↓
subduction / destruction
```

This gives the world geological memory.

At the current scale, crust age may not have a major visual effect immediately. However, it is inexpensive to add and creates a foundation for future rules.

Potential future uses include:

- ocean depth increasing with crust age
- young crust being more volcanically active
- oceanic crust becoming increasingly likely to subduct
- different erosion or terrain behavior based on geological age
- visualization/debugging of spreading patterns

---

# 16. World Data Structure

A likely initial world structure:

```rust
struct World {
    width: usize,
    height: usize,

    elevation: Vec<f32>,

    rainfall: Vec<f32>,
    water: Vec<f32>,
    flow: Vec<f32>,

    plate_id: Vec<u16>,

    crust_age: Vec<f32>,
}
```

Additional fields can be added as the simulation develops.

The large arrays should remain flat for performance and cache friendliness.

Plates themselves remain a small collection:

```rust
Vec<Plate>
```

---

# 17. Hydrology and Erosion

Erosion is a separate process from tectonics.

The intended conceptual pipeline is:

```text
Elevation
    ↓
Rainfall
    ↓
Water movement
    ↓
Flow direction
    ↓
Flow accumulation
    ↓
Erosion
    ↓
Sediment deposition
    ↓
New elevation
```

Major rivers should emerge naturally from flow accumulation rather than being explicitly generated as lines.

This is important because tectonics and erosion then interact naturally:

```text
tectonics
    ↓
mountains
    ↓
rainfall
    ↓
river systems
    ↓
erosion
    ↓
valleys / sediment
```

---

# 18. Simulation Iteration

A single simulation iteration should roughly follow this sequence:

## Step 1 — Move plates

```text
plate.position += plate.velocity
```

## Step 2 — Recalculate Voronoi diagram

Determine which plate owns each terrain point.

## Step 3 — Identify boundaries

Find neighboring plate pairs and their boundary locations.

## Step 4 — Calculate boundary normals

Determine the local direction perpendicular to each boundary.

## Step 5 — Calculate relative motion

```text
relative_velocity = velocity_a - velocity_b
```

## Step 6 — Classify boundaries

Determine:

```text
convergent
divergent
transform
```

## Step 7 — Calculate tectonic strength

Based on:

- relative velocity
- crust type
- boundary type
- local geometry

## Step 8 — Generate localized effects

Produce uplift, depression, volcanism and stress fields.

## Step 9 — Apply tectonic effects

Modify the elevation field.

## Step 10 — Apply volcanism

Generate localized volcanic features.

## Step 11 — Run hydrology/erosion

Allow the terrain to respond to rainfall and water movement.

## Step 12 — Advance geological time

Update quantities such as crust age.

---

# 19. Volcanism

Volcanoes should not simply be placed randomly along every boundary.

Instead, volcanism should be probabilistic but deterministic.

Conceptually:

```text
volcano_probability =
    base_probability
    × tectonic_strength
    × persistent_noise
```

A seeded noise/RNG system should be used so that:

```text
same seed
    → same world
```

This allows reproducible procedural generation.

Different boundary types can have different volcanic tendencies.

For example:

```text
ocean/ocean divergence
    high

continental rift
    moderate/high

ocean/continental subduction
    high

continental collision
    lower

transform
    low
```

These values are game parameters rather than strict geological probabilities.

---

# 20. Initial Terrain Rules

The original conceptual rules remain useful as the starting point.

### Convergent

Same crust:

```text
uplift both sides
```

Different crust:

```text
lower-density/appropriate side
    → depression/trench

other side
    → uplift
```

### Divergent

```text
depression around boundary
+
localized elevated volcanic/spreading features
```

### Transform

```text
small deformation
+
stress/fissures
```

These should now be implemented through localized falloff fields rather than uniformly changing every point near a plate boundary.

---

# 21. Resolution Constraint

The fixed terrain resolution is a significant design constraint.

```text
1024 × 512
```

is large enough to produce interesting continental-scale terrain, but small enough that geological features can easily become visually smeared if the deformation radius is too large.

Therefore:

> Geological effects should be narrow by default and allowed to accumulate over many iterations.

This is preferable to making each tectonic event extremely broad.

For example:

```text
iteration 1
    + small uplift

iteration 2
    + small uplift

iteration 3
    + small uplift

...

many iterations
    → major mountain chain
```

This allows geological structures to become large through **history**, rather than requiring every individual boundary calculation to modify huge regions.

---

# 22. Design Philosophy

The simulation should follow a principle of:

> **Physically motivated, but gamey.**

We are not attempting to solve plate tectonics.

Instead, each rule should have a recognizable geological motivation while remaining:

- controllable
- deterministic
- performant
- visually clear
- easy to debug
- easy to tune

The intended result is that a player looking at the generated world thinks:

> "These mountains, trenches, islands, and rivers look like they formed through geological history."

rather than:

> "This is a mathematically accurate tectonic simulation."

---

# 23. Debugging and Visualization

Because this is a simulation rather than merely a terrain generator, debugging views will be important.

Useful future debug layers include:

### Plate ownership

Each plate rendered with a distinct color.

### Plate velocity

Arrows showing `vx/vy`.

### Boundary classification

```text
convergent
divergent
transform
```

shown separately.

### Tectonic strength

Heatmap showing where deformation is occurring.

### Crust age

Heatmap from young → old crust.

### Elevation

Normal terrain visualization.

### Flow accumulation

Useful for validating river formation.

These views should make it possible to determine whether a bad-looking world is caused by:

```text
plate movement
    ↓
boundary detection
    ↓
boundary classification
    ↓
tectonic deformation
    ↓
erosion
```

rather than trying to debug everything through the final terrain alone.

---

# 24. Initial Rust Architecture

A possible initial crate organization:

```text
src/
├── lib.rs
├── world.rs
├── plate.rs
├── voronoi.rs
├── boundary.rs
├── tectonics.rs
├── hydrology.rs
├── erosion.rs
├── noise.rs
└── simulation.rs
```

This is intentionally modular.

The simulation should not become a single function containing all geological rules.

---

# 25. First Prototype

The first useful prototype does **not** need the complete system.

Recommended development sequence:

### Prototype 1 — Moving plates

Implement:

- cylindrical coordinates
- N/S walls
- Voronoi plates
- plate positions
- `vx/vy`
- movement over time

Visualization:

```text
plate boundaries
+
velocity vectors
```

### Prototype 2 — Boundary detection

Add:

- neighboring plate detection
- boundary normals
- relative velocity
- convergent/divergent/transform classification

Visualization:

```text
red    = convergent
blue   = divergent
green  = transform
```

### Prototype 3 — Tectonic deformation

Add:

- distance from boundary
- falloff
- uplift
- depression
- crust-type rules

Keep the effects deliberately small.

### Prototype 4 — Volcanism

Add:

- deterministic noise
- volcanic probability
- localized volcanic uplift

### Prototype 5 — Crust age

Add:

- creation at divergent boundaries
- aging
- eventual interaction with convergence

### Prototype 6 — Hydrology

Add:

- rainfall
- flow
- accumulation

### Prototype 7 — Erosion

Add:

- erosion
- deposition
- sediment movement

Only after these systems work independently should they be aggressively coupled.

---

# 26. Open Questions / Future Decisions

The following have intentionally **not** been finalized yet:

- exact number of plates
- initial plate distribution
- plate velocity distribution
- continental/oceanic plate ratio
- whether individual plates may contain both crust types
- exact tectonic falloff functions
- exact deformation strengths
- geological time represented by one iteration
- whether plate velocities remain constant
- how plate sites behave when reaching the N/S walls
- how oceanic crust is created and destroyed
- whether crust age belongs to points or plates
- detailed erosion algorithm
- rainfall model
- sediment transport
- whether elevation should have a sea-level/ocean representation during tectonic simulation

These should remain parameters/experiments rather than being prematurely locked down.

---

# 27. Current Decisions

For clarity, the major decisions made so far are:

| Question | Decision |
|---|---|
| World geometry | Cylinder |
| E/W | Wrap |
| N/S | Walls |
| Terrain resolution | 1024 × 512 |
| Plate representation | Moving Voronoi sites |
| Initial plate motion | Linear `vx/vy` |
| Angular velocity | Defer |
| Crust types | Oceanic / Continental |
| Boundary classification | Relative plate velocity |
| Boundary effects | Localized falloff |
| Broad deformation | Avoid |
| Tectonic effects | Intermediate effect fields |
| Volcanism | Deterministic probabilistic |
| Erosion | Separate hydrology/erosion stage |
| Crust age | Add as a lightweight geological-memory field |
| Overall goal | Physically motivated game simulation |

---

# 28. Guiding Principle

The most important constraint for the project is:

> **The geological rules need to be visible in the result.**

A technically sophisticated simulation that produces subtle or incomprehensible terrain is less useful for this game than a simplified simulation whose geological history is obvious.

The system should therefore prefer:

```text
small + local + repeated
```

over:

```text
large + global + instantaneous
```

This should allow complex geological structures to emerge from relatively simple rules while preserving the visual connection between plate boundaries and the terrain they create.