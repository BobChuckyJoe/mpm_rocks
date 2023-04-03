#[cfg(test)]
mod tests {
    use nalgebra::{Matrix3, Vector3};
    use crate::math_utils::is_point_in_triangle;
    #[test]
    fn test_point_in_triangle() {
        let a = Vector3::<f64>::new(0.0, 0.0, 0.0);
        let b = Vector3::new(1.0, 0.0, 0.0);
        let c = Vector3::new(0.0, 1.0, 0.0);
        let p = Vector3::new(0.3, 0.3, 0.0);
        assert!(is_point_in_triangle(p, a, b, c));
    }
    #[test]
    fn test_point_in_triangle_2() {
        let a = Vector3::<f64>::new(0.0, 0.0, 0.0);
        let b = Vector3::new(1.0, 0.0, 0.0);
        let c = Vector3::new(0.0, 1.0, 0.0);
        let p = Vector3::new(0.3, 1.0, 0.0);
        assert!(!is_point_in_triangle(p, a, b, c));
    }
}