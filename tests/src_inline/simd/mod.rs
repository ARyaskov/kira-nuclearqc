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
