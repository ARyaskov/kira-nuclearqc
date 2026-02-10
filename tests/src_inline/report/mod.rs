
use super::*;

#[test]
fn test_quantiles() {
    let v = vec![1.0f32, 2.0, 3.0, 4.0, 5.0];
    assert_eq!(median(&v), 3.0);
    assert_eq!(p90(&v), 5.0);
    assert_eq!(p99(&v), 5.0);
}
