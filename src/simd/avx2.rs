#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// AVX2 kernel set. Each routine splits into a SIMD head with four independent
/// f64 accumulators and a scalar tail; lane reduction is fixed-order for
/// determinism.

const LANES: usize = 32; // 4 × 8 f32, but accumulated as 4 × 4 f64

#[inline]
pub fn sum_f32_f64(values: &[f32]) -> f64 {
    let n = values.len();
    let head = (n / LANES) * LANES;
    let head_sum = if head == 0 {
        0.0
    } else {
        // SAFETY: head is a multiple of 32 ≤ n; pointer is in-bounds.
        unsafe { sum_f32_f64_avx2(values.as_ptr(), head) }
    };
    let mut tail = head_sum;
    for v in &values[head..] {
        tail += *v as f64;
    }
    tail
}

#[inline]
pub fn sum_f32(values: &[f32]) -> f32 {
    sum_f32_f64(values) as f32
}

#[inline]
pub fn max_f32(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let n = values.len();
    let head = (n / 8) * 8;
    let mut max = f32::NEG_INFINITY;
    if head != 0 {
        // SAFETY: head is a multiple of 8 ≤ n; pointer is in-bounds.
        max = unsafe { max_f32_avx2(values.as_ptr(), head) };
    }
    for v in &values[head..] {
        if *v > max {
            max = *v;
        }
    }
    if max.is_finite() { max } else { 0.0 }
}

#[inline]
pub fn entropy_f32(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let sum = sum_f32_f64(values);
    if sum <= 0.0 {
        return 0.0;
    }
    // ln() blocks autovectorisation, so the xlnx loop stays scalar.
    let inv = 1.0f64 / sum;
    let mut h = 0f64;
    for &v in values {
        if v > 0.0 {
            let p = (v as f64) * inv;
            h -= p * p.ln();
        }
    }
    h as f32
}

#[inline]
pub fn backend_name() -> &'static str {
    "avx2"
}

#[target_feature(enable = "avx2")]
unsafe fn sum_f32_f64_avx2(ptr: *const f32, len: usize) -> f64 {
    // SAFETY: caller guarantees `len` is a multiple of LANES ≤ valid length.
    unsafe {
        let mut a0 = _mm256_setzero_pd();
        let mut a1 = _mm256_setzero_pd();
        let mut a2 = _mm256_setzero_pd();
        let mut a3 = _mm256_setzero_pd();

        let mut i = 0usize;
        while i + LANES <= len {
            let v0 = _mm256_loadu_ps(ptr.add(i));
            let v1 = _mm256_loadu_ps(ptr.add(i + 8));
            let v2 = _mm256_loadu_ps(ptr.add(i + 16));
            let v3 = _mm256_loadu_ps(ptr.add(i + 24));

            let v0_lo = _mm256_cvtps_pd(_mm256_castps256_ps128(v0));
            let v0_hi = _mm256_cvtps_pd(_mm256_extractf128_ps(v0, 1));
            let v1_lo = _mm256_cvtps_pd(_mm256_castps256_ps128(v1));
            let v1_hi = _mm256_cvtps_pd(_mm256_extractf128_ps(v1, 1));
            let v2_lo = _mm256_cvtps_pd(_mm256_castps256_ps128(v2));
            let v2_hi = _mm256_cvtps_pd(_mm256_extractf128_ps(v2, 1));
            let v3_lo = _mm256_cvtps_pd(_mm256_castps256_ps128(v3));
            let v3_hi = _mm256_cvtps_pd(_mm256_extractf128_ps(v3, 1));

            a0 = _mm256_add_pd(a0, _mm256_add_pd(v0_lo, v0_hi));
            a1 = _mm256_add_pd(a1, _mm256_add_pd(v1_lo, v1_hi));
            a2 = _mm256_add_pd(a2, _mm256_add_pd(v2_lo, v2_hi));
            a3 = _mm256_add_pd(a3, _mm256_add_pd(v3_lo, v3_hi));

            i += LANES;
        }

        // Deterministic horizontal reduce: sum all 16 lanes in memory order.
        let mut buf = [0f64; 16];
        _mm256_storeu_pd(buf.as_mut_ptr(), a0);
        _mm256_storeu_pd(buf.as_mut_ptr().add(4), a1);
        _mm256_storeu_pd(buf.as_mut_ptr().add(8), a2);
        _mm256_storeu_pd(buf.as_mut_ptr().add(12), a3);
        let mut sum = 0f64;
        for &lane in &buf {
            sum += lane;
        }
        sum
    }
}

#[target_feature(enable = "avx2")]
unsafe fn max_f32_avx2(ptr: *const f32, len: usize) -> f32 {
    // SAFETY: same in-bounds guarantee as sum_f32_f64_avx2.
    unsafe {
        let mut m = _mm256_set1_ps(f32::NEG_INFINITY);
        let mut i = 0usize;
        while i + 8 <= len {
            let v = _mm256_loadu_ps(ptr.add(i));
            m = _mm256_max_ps(m, v);
            i += 8;
        }
        let mut buf = [0f32; 8];
        _mm256_storeu_ps(buf.as_mut_ptr(), m);
        let mut out = f32::NEG_INFINITY;
        for &lane in &buf {
            if lane > out {
                out = lane;
            }
        }
        out
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/simd/avx2.rs"]
mod tests;
