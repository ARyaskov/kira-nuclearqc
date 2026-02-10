#[inline]
pub fn sum_f32_f64(values: &[f32]) -> f64 {
    backend::sum_f32_f64(values)
}

#[inline]
pub fn sum_f32(values: &[f32]) -> f32 {
    sum_f32_f64(values) as f32
}

#[inline]
pub fn max_f32(values: &[f32]) -> f32 {
    backend::max_f32(values)
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
    let mut h = 0f64;
    for &v in values {
        let p = (v as f64) / sum;
        if p > 0.0 {
            h -= p * p.ln();
        }
    }
    h as f32
}

#[inline]
pub fn backend_name() -> &'static str {
    backend::backend_name()
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
mod backend {
    pub use crate::simd::avx2::*;
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
mod backend {
    pub use crate::simd::neon::*;
}

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    all(target_arch = "aarch64", target_feature = "neon"),
)))]
mod backend {
    pub use crate::simd::scalar::*;
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub mod avx2;
#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
pub mod neon;
pub mod scalar;

#[cfg(test)]
#[path = "../../tests/src_inline/simd/mod.rs"]
mod tests;
