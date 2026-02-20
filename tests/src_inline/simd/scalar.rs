use super::*;

#[test]
fn test_sum() {
    let v = [1.0f32, 2.0, 3.0];
    assert_eq!(sum_f32_f64(&v), 6.0);
}

#[test]
fn test_max() {
    let v = [1.0f32, -1.0, 3.0];
    assert_eq!(max_f32(&v), 3.0);
}
