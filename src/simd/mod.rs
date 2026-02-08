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
mod tests {
    use super::*;

    #[test]
    fn test_backend_name() {
        let name = backend_name();
        #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
        assert_eq!(name, "avx2");
        #[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
        assert_eq!(name, "neon");
        #[cfg(not(any(
            all(target_arch = "x86_64", target_feature = "avx2"),
            all(target_arch = "aarch64", target_feature = "neon"),
        )))]
        assert_eq!(name, "scalar");
    }

    #[test]
    fn test_sum_determinism() {
        let values = [0.1f32, 0.2, 0.3, 0.4, 0.5, 0.6];
        let a = sum_f32_f64(&values);
        let b = sum_f32_f64(&values);
        assert!(a.to_bits() == b.to_bits());
    }

    #[test]
    fn test_backend_equiv_scalar() {
        let values = [0.1f32, 0.2, 0.3, 0.4, 0.5, 0.6, 1.5];
        assert_eq!(sum_f32_f64(&values), scalar::sum_f32_f64(&values));
        assert_eq!(max_f32(&values), scalar::max_f32(&values));
        assert_eq!(entropy_f32(&values), scalar::entropy_f32(&values));
    }
}
