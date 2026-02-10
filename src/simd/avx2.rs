#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

pub fn sum_f32_f64(values: &[f32]) -> f64 {
    // Deterministic order: process chunks, but accumulate each lane in order.
    let mut sum = 0f64;
    let mut i = 0usize;
    let n = values.len();
    unsafe {
        while i + 8 <= n {
            let ptr = values.as_ptr().add(i) as *const __m256;
            let v = _mm256_loadu_ps(ptr);
            let mut lanes = [0f32; 8];
            _mm256_storeu_ps(lanes.as_mut_ptr(), v);
            for lane in &lanes {
                sum += *lane as f64;
            }
            i += 8;
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
        while i + 8 <= n {
            let ptr = values.as_ptr().add(i) as *const __m256;
            let v = _mm256_loadu_ps(ptr);
            let mut lanes = [0f32; 8];
            _mm256_storeu_ps(lanes.as_mut_ptr(), v);
            for lane in &lanes {
                if *lane > max {
                    max = *lane;
                }
            }
            i += 8;
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
    "avx2"
}

#[cfg(test)]
#[path = "../../tests/src_inline/simd/avx2.rs"]
mod tests;
