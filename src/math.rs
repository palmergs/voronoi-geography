//! Cylinder math helpers.
//!
//! Vector types come from `glam`; this module only adds the wrapping rules the
//! cylindrical world needs (x wraps, y is walled).

pub use glam::Vec2;

/// Shortest signed distance from `b` to `a` on a periodic axis of `period`.
///
/// The result is always in `[-period/2, period/2]`.
#[inline]
pub fn wrap_delta(a: f32, b: f32, period: f32) -> f32 {
    let half = period * 0.5;
    let d = (a - b).rem_euclid(period);
    if d > half { d - period } else { d }
}

/// 90 degree rotation, used for boundary tangents.
#[inline]
pub fn perp(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}

#[inline]
pub fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

/// Hermite ramp from 0 at `edge0` to 1 at `edge1`.
#[inline]
pub fn smoothstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    if (edge1 - edge0).abs() < 1e-9 {
        return if x < edge0 { 0.0 } else { 1.0 };
    }
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// A gaussian bump of height 1 centred on `center`.
#[inline]
pub fn bump(x: f32, center: f32, sigma: f32) -> f32 {
    let t = (x - center) / sigma.max(1e-6);
    (-0.5 * t * t).exp()
}

/// Value at a given quantile of `data`.
///
/// Used to auto-scale things to a world's own range without letting a single
/// outlier cell flatten everything else - a river threshold, a colour ramp.
/// Runs in O(n) and ignores non-finite values.
pub fn percentile(data: &[f32], q: f32) -> f32 {
    if data.is_empty() {
        return 0.0;
    }
    let mut v: Vec<f32> = data.iter().copied().filter(|f| f.is_finite()).collect();
    if v.is_empty() {
        return 0.0;
    }
    let k = ((v.len() - 1) as f32 * q.clamp(0.0, 1.0)).round() as usize;
    let (_, nth, _) = v.select_nth_unstable_by(k, f32::total_cmp);
    *nth
}

/// Unit vector, or zero for a degenerate input.
#[inline]
pub fn normalize_or_zero(v: Vec2) -> Vec2 {
    v.normalize_or_zero()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrap_delta_takes_the_short_way_around() {
        assert_eq!(wrap_delta(1.0, 1023.0, 1024.0), 2.0);
        assert_eq!(wrap_delta(1023.0, 1.0, 1024.0), -2.0);
        assert_eq!(wrap_delta(10.0, 4.0, 1024.0), 6.0);
    }

    #[test]
    fn percentile_ignores_outliers() {
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();
        assert_eq!(percentile(&data, 0.0), 0.0);
        assert_eq!(percentile(&data, 0.5), 50.0);
        assert_eq!(percentile(&data, 1.0), 99.0);
        assert_eq!(percentile(&[], 0.5), 0.0);
        assert_eq!(percentile(&[f32::NAN, 3.0], 1.0), 3.0);
    }

    #[test]
    fn smoothstep_and_bump_behave() {
        assert_eq!(smoothstep(0.0, 1.0, -1.0), 0.0);
        assert_eq!(smoothstep(0.0, 1.0, 2.0), 1.0);
        assert!((smoothstep(0.0, 1.0, 0.5) - 0.5).abs() < 1e-6);
        assert!((bump(3.0, 3.0, 2.0) - 1.0).abs() < 1e-6);
        assert!(bump(9.0, 3.0, 1.0) < 0.01);
    }

    #[test]
    fn wrap_delta_stays_in_half_period() {
        for a in 0..1024 {
            let d = wrap_delta(a as f32, 0.0, 1024.0);
            assert!((-512.0..=512.0).contains(&d), "a={a} d={d}");
        }
    }
}
