pub fn sum_f32_f64(values: &[f32]) -> f64 {
    let mut sum = 0f64;
    for &v in values {
        sum += v as f64;
    }
    sum
}

pub fn max_f32(values: &[f32]) -> f32 {
    let mut max = f32::NEG_INFINITY;
    for &v in values {
        if v > max {
            max = v;
        }
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
    "scalar"
}

#[cfg(test)]
mod tests {
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
}
