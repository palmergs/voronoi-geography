# voronoi-geography

A game-oriented plate tectonics world simulation, implementing
[`docs/plate-tectonics-world-sim-project-design.md`](docs/plate-tectonics-world-sim-project-design.md).

The world is a `1024 x 512` cylindrical heightfield: east/west wraps, north/south
are walls. Tectonic plates are moving Voronoi sites that drift across the fixed
grid. Boundary behaviour is derived from relative plate velocity, and every
effect is applied through a narrow distance falloff so features stay legible at
this resolution. Geology accumulates over hundreds of small steps rather than a
few large ones.

## Quick start

```sh
cargo run --release -- --seed 7 --steps 400 --layers all --out out
```

```
--seed <N>           World seed (default 1)
--steps <N>          Simulation iterations (default 240)
--plates <N>         Number of tectonic plates (default 32)
--continental <F>    Fraction of the world that is continental crust (default 0.4)
--width/--height     World size (default 1024 x 512)
--dt <F>             Geological time per iteration (default 1.0)
--erosion-every <N>  Run hydrology/erosion every N steps (default 2, 0 = off)
--layers <A,B,...>   Layers to write, or 'all'
--snapshot-every <N> Also write numbered snapshots, for watching the history
```

Roughly 57 ms/step at full resolution with 32 plates, so a 400-step world takes
about 25 seconds.

## Layers

Each layer isolates one stage, so a bad-looking world can be traced to the stage
that caused it rather than debugged through the final terrain.

| Layer | Shows |
|---|---|
| `terrain` | Hypsometric tint, hillshade, rivers and lakes |
| `elevation` | Elevation alone, unshaded |
| `plates` | Plate ownership, sites, velocity arrows |
| `boundaries` | Classification: red convergent, blue divergent, green transform |
| `strength` | Where deformation is being applied, and how hard |
| `crust-age` | Young to old crust - shows seafloor spreading directly |
| `rainfall` | Latitude bands, orographic gain and rain shadow |
| `flow` | Flow accumulation on a log scale |

## Modules

| Module | Role | Design doc |
|---|---|---|
| `world` | The fixed grid and its cylindrical index rules | 2 |
| `plate` | Plates, the crust field, plate motion | 3, 4, 14 |
| `voronoi` | Ownership each step, plus the boundary warp | 5 |
| `boundary` | Detection, classification, distance field | 6-9 |
| `tectonics` | Deformation rules, falloff, volcanism | 7-14, 19-21 |
| `hydrology` | Rainfall, depression filling, flow routing | 17 |
| `erosion` | Incision, sediment transport, hillslope creep | 17 |
| `simulation` | The iteration order and the parameter set | 18 |
| `render` | Debug and terrain layers | 23 |
| `noise`, `math` | Seeded noise, cylinder math | - |

Everything is deterministic: the same seed always produces the same world.

## Decisions

The design doc (section 26) deliberately left a number of questions open. What
this implementation settled on, and why:

- **Plates at the walls reflect.** A site that would walk off the top or bottom
  bounces. Clamping would pile plates up against the wall; wrapping would
  contradict the cylinder.
- **Crust type comes from a noise field, not from plate ownership.** If crust
  type were read off the plate, every coastline would be a Voronoi edge and the
  world would look like polygon soup. A large-scale `CrustField` decides what is
  continent; plates take their type from the field under their site, and may
  straddle both kinds of crust.
- **Crust age lives on terrain cells.** Divergent boundaries reset it to zero and
  it grows everywhere else, so the crust-age layer reads as a spreading pattern.
  The per-plate mean is used for one decision: which of two oceanic plates
  subducts (the older, denser one).
- **Boundary effects use nearest-boundary attribution.** Each cell is affected by
  the single closest boundary, not the sum of every boundary in range. This is
  what keeps a mountain range traceable to one boundary.
- **Falloff is `exp(-d / 0.35r)`**, tapered over the last 30% of the radius and
  hard-zero beyond `r`, so the cutoff never shows up as a visible circular edge.
- **Sea level is fixed at 0.0** for the whole run.
- **Plate velocities are constant**, apart from wall reflection. Linear `vx/vy`
  only; no Euler poles.

### One addition beyond the design

Plate boundaries are **domain-warped** before the Voronoi lookup
(`voronoi::Warp`). Straight perpendicular bisectors produced dead-straight
mountain belts and left the polygon pattern visible through the finished world.
The warp is a fixed, precomputed distortion of the lookup position, so the
diagram stays deterministic and seamless while boundaries bend into something
that reads as a rift or a suture.

### Notes

- `CylinderNoise` maps x onto a circle in 3D so the noise is seamless across the
  x=0 seam by construction, rather than by a periodicity hack.
- Noise fields are **calibrated at construction** to a known standard deviation.
  Perlin fBm only spans about ±0.6 in practice, and by a factor that depends on
  the octave count, so without this an amplitude parameter would not mean what
  it says.
- The boundary distance field is an exact Euclidean distance transform
  (Felzenszwalb & Huttenlocher), O(cells), with the wrap handled by running the
  row pass over three tiled copies. A chamfer approximation was tried first and
  its ~5% error near the seam was too coarse for the falloffs to key off.
- Flow routing is priority-flood with an epsilon tilt, which fills every
  depression and yields a drainage direction for every cell in one pass. The
  order cells are popped in is a topological order of the network, so flow
  accumulation and sediment transport are single reverse sweeps.

## Tests

```sh
cargo test
```

76 tests. The interesting ones assert geological behaviour rather than
implementation details: that continental collision lifts both sides while
ocean/continent convergence digs a trench on one side and builds a cordillera on
the other, that a mid-ocean ridge stands above its flanks while a continental
rift sits below its shoulders, that effects never reach plate interiors, that
water is conserved through the drainage network, and that the same seed gives a
byte-identical world.
