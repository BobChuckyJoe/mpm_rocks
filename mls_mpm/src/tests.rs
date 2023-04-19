#[cfg(test)]
mod tests {
    use nalgebra::{Matrix3, Vector3};
    use crate::config::{GRID_SPACING, GRID_LENGTH};
    use crate::{math_utils::is_point_in_triangle, config::SIMULATION_SIZE};
    use crate::equations::{get_base_grid_ind, weighting_function};
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
    
    #[test]
    fn test_weighting_function_sums_to_unity() {
        let particle_position = Vector3::<f64>::new(1.0, 1.0, 1.0) * SIMULATION_SIZE / 2.0;
        let particle_base_coord = get_base_grid_ind(&particle_position, GRID_SPACING);  
        
        let mut weighting_sum = 0.0;
        for dx in -2..3 {
            for dy in -2..3 {
                for dz in -2..3 {
                    let x = particle_base_coord.0 as i64 + dx;
                    let y = particle_base_coord.1 as i64 + dy;
                    let z = particle_base_coord.2 as i64 + dz;
                    if x < 0
                        || x >= GRID_LENGTH as i64
                        || y < 0
                        || y >= GRID_LENGTH as i64
                        || z < 0
                        || z >= GRID_LENGTH as i64
                    {
                        continue;
                    }
        
                    let x = x as usize;
                    let y = y as usize;
                    let z = z as usize;
                    weighting_sum += weighting_function(particle_position, (x, y, z));
                }
            }
        }
        assert!((weighting_sum - 1.0).abs() < 1e-6, "weighting sum is {}", weighting_sum);
    }
}