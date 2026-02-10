
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
