//! Command line driver: run the simulation and write out debug layers.

use std::path::{Path, PathBuf};
use std::time::Instant;

use voronoi_geology::plate::PlateParams;
use voronoi_geology::render::{self, Layer};
use voronoi_geology::simulation::{Simulation, SimulationParams};
use voronoi_geology::world::{DEFAULT_HEIGHT, DEFAULT_WIDTH};

struct Args {
    params: SimulationParams,
    steps: u64,
    out: PathBuf,
    layers: Vec<Layer>,
    /// Also write a numbered snapshot every N steps, for watching the history.
    snapshot_every: u64,
    quiet: bool,
}

impl Default for Args {
    fn default() -> Self {
        Args {
            params: SimulationParams::default(),
            steps: 240,
            out: PathBuf::from("out"),
            layers: vec![
                Layer::Terrain,
                Layer::Plates,
                Layer::Boundaries,
                Layer::CrustAge,
                Layer::Flow,
            ],
            snapshot_every: 0,
            quiet: false,
        }
    }
}

const USAGE: &str = "\
voronoi-geology - plate tectonics world simulation

USAGE:
    voronoi-geology [OPTIONS]

OPTIONS:
    --seed <N>           World seed (default clock millis)
    --steps <N>          Simulation iterations (default 240)
    --plates <N>         Number of tectonic plates (default 32)
    --continental <F>    Fraction of plates that are continental (default 0.4)
    --width <N>          World width, wraps east/west (default 1024)
    --height <N>         World height, walled north/south (default 512)
    --dt <F>             Geological time per iteration (default 1.0)
    --sea-level <F>      Where the sea sits (default 0.0). Lower it to drain
                         the world without changing the tectonics
    --erosion-every <N>  Run hydrology/erosion every N steps (default 2, 0 = off)
    --out <DIR>          Output directory (default ./out)
    --layers <A,B,...>   Layers to write, or 'all' (default terrain,plates,
                         boundaries,crust-age,flow)
    --snapshot-every <N> Also write numbered snapshots every N steps
    --quiet              Only print the final summary
    --help               Show this message

LAYERS:
    terrain elevation plates boundaries strength crust-age rainfall flow
";

fn main() {
    let args = match parse_args() {
        Ok(Some(args)) => args,
        Ok(None) => return,
        Err(e) => {
            eprintln!("error: {e}\n\n{USAGE}");
            std::process::exit(2);
        }
    };

    if let Err(e) = std::fs::create_dir_all(&args.out) {
        eprintln!("error: cannot create {}: {e}", args.out.display());
        std::process::exit(1);
    }

    let p = &args.params;
    if !args.quiet {
        println!(
            "seed {} | {}x{} | {} plates ({:.0}% continental) | {} steps",
            p.seed,
            p.width,
            p.height,
            p.plates.count,
            p.plates.continental_fraction * 100.0,
            args.steps
        );
    }

    let started = Instant::now();
    let mut sim = Simulation::new(args.params);
    if !args.quiet {
        println!("  setup {:>6.2}s", started.elapsed().as_secs_f32());
    }

    let sim_start = Instant::now();
    for step in 1..=args.steps {
        sim.step();

        if !args.quiet && (step % 20 == 0 || step == args.steps) {
            let s = &sim.stats;
            println!(
                "  step {step:>4}  conv {:>5} div {:>5} trans {:>5}  volc {:>3}  \
                 elev {:>6.2}..{:<5.2}  land {:>4.1}%",
                s.convergent,
                s.divergent,
                s.transform,
                s.volcanoes,
                s.min_elevation,
                s.max_elevation,
                s.land_fraction * 100.0
            );
        }

        if args.snapshot_every > 0 && step % args.snapshot_every == 0 {
            write_layers(&sim, &args.layers, &args.out, Some(step), args.quiet);
        }
    }
    let sim_time = sim_start.elapsed();

    write_layers(&sim, &args.layers, &args.out, None, args.quiet);

    let s = &sim.stats;
    println!(
        "\n{} steps in {:.2}s ({:.0} ms/step)",
        args.steps,
        sim_time.as_secs_f32(),
        sim_time.as_secs_f32() * 1000.0 / args.steps.max(1) as f32
    );
    println!(
        "elevation {:.2} to {:.2} | land {:.1}% | boundaries {} conv / {} div / {} trans",
        s.min_elevation,
        s.max_elevation,
        s.land_fraction * 100.0,
        s.convergent,
        s.divergent,
        s.transform
    );
    println!(
        "erosion moved {:.0} units, deposited {:.0}, exported {:.0}",
        sim.totals.eroded, sim.totals.deposited, sim.totals.exported
    );
    println!("wrote {} layers to {}", args.layers.len(), args.out.display());
}

fn write_layers(sim: &Simulation, layers: &[Layer], out: &Path, step: Option<u64>, quiet: bool) {
    for layer in layers {
        let name = match step {
            Some(n) => format!("{}-{:05}.png", layer.name(), n),
            None => format!("{}.png", layer.name()),
        };
        let path = out.join(&name);
        let img = render::render(sim, *layer);
        if let Err(e) = img.save(&path) {
            eprintln!("error: cannot write {}: {e}", path.display());
        } else if !quiet && step.is_none() {
            println!("  wrote {}", path.display());
        }
    }
}

fn parse_args() -> Result<Option<Args>, String> {
    let mut args = Args::default();
    let mut argv = std::env::args().skip(1);

    while let Some(arg) = argv.next() {
        let mut value = || -> Result<String, String> {
            argv.next()
                .ok_or_else(|| format!("{arg} needs a value"))
        };

        match arg.as_str() {
            "--help" | "-h" => {
                print!("{USAGE}");
                return Ok(None);
            }
            "--quiet" | "-q" => args.quiet = true,
            "--seed" => args.params.seed = parse(&value()?, &arg)?,
            "--steps" => args.steps = parse(&value()?, &arg)?,
            "--plates" => args.params.plates.count = parse(&value()?, &arg)?,
            "--continental" => {
                args.params.plates.continental_fraction = parse(&value()?, &arg)?
            }
            "--width" => args.params.width = parse(&value()?, &arg)?,
            "--height" => args.params.height = parse(&value()?, &arg)?,
            "--dt" => args.params.dt = parse(&value()?, &arg)?,
            "--sea-level" => args.params.sea_level = parse(&value()?, &arg)?,
            "--erosion-every" => args.params.erosion_interval = parse(&value()?, &arg)?,
            "--out" => args.out = PathBuf::from(value()?),
            "--snapshot-every" => args.snapshot_every = parse(&value()?, &arg)?,
            "--layers" => {
                let raw = value()?;
                args.layers = if raw == "all" {
                    Layer::ALL.to_vec()
                } else {
                    raw.split(',')
                        .map(|n| {
                            Layer::parse(n.trim())
                                .ok_or_else(|| format!("unknown layer '{}'", n.trim()))
                        })
                        .collect::<Result<_, _>>()?
                };
            }
            other => return Err(format!("unknown option '{other}'")),
        }
    }

    validate(&args)?;
    Ok(Some(args))
}

fn parse<T: std::str::FromStr>(raw: &str, flag: &str) -> Result<T, String> {
    raw.parse()
        .map_err(|_| format!("{flag}: '{raw}' is not a valid value"))
}

fn validate(args: &Args) -> Result<(), String> {
    let p = &args.params;
    if p.width == 0 || p.height < 3 {
        return Err("world must be at least 1x3".into());
    }
    if p.plates.count < 2 {
        return Err("need at least 2 plates for there to be any boundaries".into());
    }
    if p.plates.count > u16::MAX as usize - 1 {
        return Err(format!("at most {} plates", u16::MAX - 1));
    }
    if !(0.0..=1.0).contains(&p.plates.continental_fraction) {
        return Err("--continental must be between 0 and 1".into());
    }
    if p.dt <= 0.0 {
        return Err("--dt must be positive".into());
    }
    if args.layers.is_empty() {
        return Err("no layers selected".into());
    }
    let _ = (DEFAULT_WIDTH, DEFAULT_HEIGHT, PlateParams::default());
    Ok(())
}
