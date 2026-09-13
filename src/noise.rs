//! Deterministic noise and randomness.
//!
//! Everything must be reproducible from a seed, so there are no global RNGs.
//!
//! * [`rng`] returns a seeded `rand_pcg` stream for setup work (plate placement).
//! * [`CylinderNoise`] wraps the `noise` crate's fBm. Because the world is a
//!   cylinder, x is mapped onto a circle in 3D and the noise is sampled on that
//!   circle - the field is then seamless across the x=0 seam by construction,
//!   with no periodicity hacks.
//! * [`hash_unit`] is a stateless coordinate hash, for per-cell decisions that
//!   must not depend on iteration order (volcano placement, transform jitter).

use noise::{Fbm, MultiFractal, NoiseFn, Perlin};
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;

pub type Rng = Pcg64Mcg;

/// A reproducible RNG stream. `stream` separates independent uses of one seed.
pub fn rng(seed: u64, stream: u64) -> Rng {
    Pcg64Mcg::seed_from_u64(seed ^ stream.wrapping_mul(0x9E37_79B9_7F4A_7C15))
}

/// Seamless fractal noise over the cylindrical world.
pub struct CylinderNoise {
    fbm: Fbm<Perlin>,
    /// Radius of the circle x is mapped onto, in noise units.
    radius: f64,
    /// World width, for the x -> angle mapping.
    width: f64,
    /// Scale applied to y so features stay isotropic.
    y_scale: f64,
    /// Calibration factor, see [`CylinderNoise::new`].
    gain: f32,
}

/// Target standard deviation of a calibrated field. Chosen so that samples
/// mostly land in [-1, 1] with occasional excursions past it, which is what
/// callers assume when they write an amplitude in world units.
const TARGET_SD: f32 = 0.32;

impl CylinderNoise {
    /// `feature_size` is roughly the width of one noise feature, in cells.
    pub fn new(seed: u32, width: usize, feature_size: f32, octaves: usize) -> Self {
        let fbm = Fbm::<Perlin>::new(seed)
            .set_octaves(octaves)
            .set_persistence(0.5)
            .set_lacunarity(2.0);
        // One full trip around the cylinder must cover width/feature_size
        // noise units, so the circle's circumference is exactly that.
        let circumference = width as f64 / feature_size as f64;
        let mut field = Self {
            fbm,
            radius: circumference / std::f64::consts::TAU,
            width: width as f64,
            y_scale: 1.0 / feature_size as f64,
            gain: 1.0,
        };

        // Perlin fBm does not fill [-1, 1]: with these octave counts it comes
        // out around a third of that, and the amount depends on how many
        // octaves were asked for. Left uncalibrated, every amplitude parameter
        // in the simulation would quietly mean something different from what
        // it says. So measure the field once and scale it to a known spread.
        let mut sum = 0.0f64;
        let mut sum_sq = 0.0f64;
        let samples = 2048;
        for i in 0..samples {
            // A low-discrepancy walk over the field, so the estimate does not
            // depend on lining up with the noise lattice.
            let t = i as f32;
            let x = (t * 0.754_877_7).fract() * width as f32;
            let y = (t * 0.569_840_3).fract() * (feature_size * 64.0);
            let v = field.raw(x, y);
            sum += v as f64;
            sum_sq += (v * v) as f64;
        }
        let mean = sum / samples as f64;
        let sd = (sum_sq / samples as f64 - mean * mean).max(1e-12).sqrt() as f32;
        field.gain = (TARGET_SD / sd).clamp(0.5, 8.0);
        field
    }

    /// The underlying field, before calibration.
    fn raw(&self, x: f32, y: f32) -> f32 {
        let theta = (x as f64 / self.width) * std::f64::consts::TAU;
        self.fbm.get([
            self.radius * theta.cos(),
            self.radius * theta.sin(),
            y as f64 * self.y_scale,
        ]) as f32
    }

    /// Sample at a world cell. Calibrated so the output is roughly `[-1, 1]`.
    pub fn get(&self, x: f32, y: f32) -> f32 {
        self.raw(x, y) * self.gain
    }
}

/// Stateless hash of two coordinates (plus a salt) to `[0, 1)`.
///
/// splitmix64 finaliser - cheap enough to call per cell per iteration.
pub fn hash_unit(seed: u64, a: i64, b: i64) -> f32 {
    let mut z = seed
        .wrapping_mul(0xD6E8_FEB8_6659_FD93)
        .wrapping_add((a as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
        .wrapping_add((b as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F));
    z ^= z >> 30;
    z = z.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z ^= z >> 27;
    z = z.wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    ((z >> 40) as f32) / (1u32 << 24) as f32
}

/// Stateless coordinate hash to `[-1, 1)`.
pub fn hash_signed(seed: u64, a: i64, b: i64) -> f32 {
    hash_unit(seed, a, b) * 2.0 - 1.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::RngExt as _;

    #[test]
    fn rng_is_deterministic_and_stream_separated() {
        let mut a = rng(42, 0);
        let mut b = rng(42, 0);
        let mut c = rng(42, 1);
        let xs: Vec<u64> = (0..16).map(|_| a.random()).collect();
        let ys: Vec<u64> = (0..16).map(|_| b.random()).collect();
        let zs: Vec<u64> = (0..16).map(|_| c.random()).collect();
        assert_eq!(xs, ys);
        assert_ne!(xs, zs);
    }

    #[test]
    fn cylinder_noise_is_seamless() {
        let n = CylinderNoise::new(7, 1024, 64.0, 4);
        for i in 0..32 {
            let y = i as f32 * 13.0;
            let a = n.get(0.0, y);
            let b = n.get(1024.0, y);
            assert!((a - b).abs() < 1e-5, "seam mismatch at y={y}: {a} vs {b}");
            // Adjacent samples across the seam should also be continuous.
            let l = n.get(1023.0, y);
            assert!((l - a).abs() < 0.2, "discontinuity at seam: {l} vs {a}");
        }
    }

    #[test]
    fn cylinder_noise_actually_varies() {
        let n = CylinderNoise::new(7, 1024, 64.0, 4);
        let samples: Vec<f32> = (0..64).map(|i| n.get(i as f32 * 16.0, 128.0)).collect();
        let min = samples.iter().cloned().fold(f32::MAX, f32::min);
        let max = samples.iter().cloned().fold(f32::MIN, f32::max);
        assert!(max - min > 0.2, "noise is too flat: {min}..{max}");
    }

    #[test]
    fn fields_are_calibrated_to_a_known_spread() {
        // Whatever the feature size or octave count, an amplitude of 1.0 in
        // world units should mean about 1.0 of terrain.
        for (feature, octaves) in [(12.0f32, 3usize), (55.0, 5), (260.0, 3)] {
            let n = CylinderNoise::new(11, 1024, feature, octaves);
            let mut values = Vec::new();
            for y in (0..512).step_by(7) {
                for x in (0..1024).step_by(7) {
                    values.push(n.get(x as f32, y as f32));
                }
            }
            let mean = values.iter().sum::<f32>() / values.len() as f32;
            let sd = (values.iter().map(|v| (v - mean).powi(2)).sum::<f32>()
                / values.len() as f32)
                .sqrt();
            let peak = values.iter().fold(0.0f32, |a, v| a.max(v.abs()));

            assert!(
                (sd - TARGET_SD).abs() < 0.08,
                "feature {feature}/{octaves} octaves: sd {sd} is off target"
            );
            assert!(
                (0.7..2.0).contains(&peak),
                "feature {feature}/{octaves} octaves: peak {peak} is implausible"
            );
        }
    }

    #[test]
    fn hash_is_stable_and_in_range() {
        for a in 0..50 {
            for b in 0..50 {
                let v = hash_unit(99, a, b);
                assert!((0.0..1.0).contains(&v));
                assert_eq!(v, hash_unit(99, a, b));
            }
        }
        assert_ne!(hash_unit(1, 2, 3), hash_unit(2, 2, 3));
    }
}
