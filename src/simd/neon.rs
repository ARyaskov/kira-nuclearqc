#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

pub fn sum_f32_f64(values: &[f32]) -> f64 {
    let mut sum = 0f64;
    let mut i = 0usize;
    let n = values.len();
    unsafe {
        while i + 4 <= n {
            let ptr = values.as_ptr().add(i);
            let v = vld1q_f32(ptr);
            let mut lanes = [0f32; 4];
            vst1q_f32(lanes.as_mut_ptr(), v);
            for lane in &lanes {
                sum += *lane as f64;
            }
            i += 4;
        }
    }
    while i < n {
        sum += values[i] as f64;
        i += 1;
    }
    sum
}

pub fn max_f32(values: &[f32]) -> f32 {
    let mut max = f32::NEG_INFINITY;
    let mut i = 0usize;
    let n = values.len();
    unsafe {
        while i + 4 <= n {
            let ptr = values.as_ptr().add(i);
            let v = vld1q_f32(ptr);
            let mut lanes = [0f32; 4];
            vst1q_f32(lanes.as_mut_ptr(), v);
            for lane in &lanes {
                if *lane > max {
                    max = *lane;
                }
            }
            i += 4;
        }
    }
    while i < n {
        let v = values[i];
        if v > max {
            max = v;
        }
        i += 1;
    }
    if max.is_finite() { max } else { 0.0 }
}

pub fn entropy_f32(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let sum = sum_f32_f64(values);
    if sum <= 0.0 {
        return 0.0;
    }
    let mut h = 0f64;
    for &v in values {
        let p = (v as f64) / sum;
        if p > 0.0 {
            h -= p * p.ln();
        }
    }
    h as f32
}

pub fn backend_name() -> &'static str {
    "neon"
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simd::scalar;

    #[test]
    fn test_sum_equiv() {
        let v = [0.1f32, 0.2, 0.3, 0.4, 0.5, 0.6, 1.1];
        assert_eq!(sum_f32_f64(&v), scalar::sum_f32_f64(&v));
    }

    #[test]
    fn test_max_equiv() {
        let v = [0.1f32, 2.0, 0.3, -4.0];
        assert_eq!(max_f32(&v), scalar::max_f32(&v));
    }

    #[test]
    fn test_entropy_equiv() {
        let v = [0.1f32, 0.2, 0.3, 0.4];
        assert_eq!(entropy_f32(&v), scalar::entropy_f32(&v));
    }
}
