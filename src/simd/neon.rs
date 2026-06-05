#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

/// NEON kernel set. Mirrors avx2.rs: a SIMD head with four independent f64
/// accumulators and a scalar tail, fixed-order lane reduction.

const LANES: usize = 16; // 4 × 4 f32, accumulated as 4 × 2 f64

#[inline]
pub fn sum_f32_f64(values: &[f32]) -> f64 {
    let n = values.len();
    let head = (n / LANES) * LANES;
    let head_sum = if head == 0 {
        0.0
    } else {
        // SAFETY: head ≤ n; pointer is in-bounds.
        unsafe { sum_f32_f64_neon(values.as_ptr(), head) }
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
    let head = (n / 4) * 4;
    let mut max = f32::NEG_INFINITY;
    if head != 0 {
        // SAFETY: head is a multiple of 4 ≤ n; pointer is in-bounds.
        max = unsafe { max_f32_neon(values.as_ptr(), head) };
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
    "neon"
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn sum_f32_f64_neon(ptr: *const f32, len: usize) -> f64 {
    // SAFETY: caller guarantees `len ≤ ptr-valid bytes / 4`.
    unsafe {
        let mut a0 = vdupq_n_f64(0.0);
        let mut a1 = vdupq_n_f64(0.0);
        let mut a2 = vdupq_n_f64(0.0);
        let mut a3 = vdupq_n_f64(0.0);

        let mut i = 0usize;
        while i + LANES <= len {
            let v0 = vld1q_f32(ptr.add(i));
            let v1 = vld1q_f32(ptr.add(i + 4));
            let v2 = vld1q_f32(ptr.add(i + 8));
            let v3 = vld1q_f32(ptr.add(i + 12));

            a0 = vaddq_f64(a0, vcvt_f64_f32(vget_low_f32(v0)));
            a0 = vaddq_f64(a0, vcvt_high_f64_f32(v0));
            a1 = vaddq_f64(a1, vcvt_f64_f32(vget_low_f32(v1)));
            a1 = vaddq_f64(a1, vcvt_high_f64_f32(v1));
            a2 = vaddq_f64(a2, vcvt_f64_f32(vget_low_f32(v2)));
            a2 = vaddq_f64(a2, vcvt_high_f64_f32(v2));
            a3 = vaddq_f64(a3, vcvt_f64_f32(vget_low_f32(v3)));
            a3 = vaddq_f64(a3, vcvt_high_f64_f32(v3));

            i += LANES;
        }

        let mut buf = [0f64; 8];
        vst1q_f64(buf.as_mut_ptr(), a0);
        vst1q_f64(buf.as_mut_ptr().add(2), a1);
        vst1q_f64(buf.as_mut_ptr().add(4), a2);
        vst1q_f64(buf.as_mut_ptr().add(6), a3);
        let mut sum = 0f64;
        for &lane in &buf {
            sum += lane;
        }
        sum
    }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn max_f32_neon(ptr: *const f32, len: usize) -> f32 {
    // SAFETY: same in-bounds guarantee as sum_f32_f64_neon.
    unsafe {
        let mut m = vdupq_n_f32(f32::NEG_INFINITY);
        let mut i = 0usize;
        while i + 4 <= len {
            let v = vld1q_f32(ptr.add(i));
            m = vmaxq_f32(m, v);
            i += 4;
        }
        let mut buf = [0f32; 4];
        vst1q_f32(buf.as_mut_ptr(), m);
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
#[path = "../../tests/src_inline/simd/neon.rs"]
mod tests;
