/// Random helper functions for math things
use nalgebra::{Matrix3, Vector3, UnitQuaternion};

// Helper function to iterate through 3x3x3 grid around the starting index
// This is mainly used by rigid particle splatting
pub fn iterate_over_3x3(start: (usize, usize, usize), grid_lengths: (usize, usize, usize)) -> Vec<(usize, usize, usize)> {
    let mut inds = Vec::new();
    let (i, j, k) = start;
    for di in -1..2 {
        for dj in -1..2 {
            for dk in -1..2 {
                let new_i = i as i32 + di;
                let new_j = j as i32 + dj;
                let new_k = k as i32 + dk;
                if (new_i >= 0 && new_i < grid_lengths.0 as i32)
                    && (new_j >= 0 && new_j < grid_lengths.1 as i32)
                    && (new_k >= 0 && new_k < grid_lengths.2 as i32)
                {
                    inds.push((new_i as usize, new_j as usize, new_k as usize));
                }
            }
        }
    }
    assert!(inds.len() <= 27);
    inds
}

pub fn project_point_into_plane(
    point: Vector3<f64>,
    plane_normal: Vector3<f64>,
    plane_point: Vector3<f64>,
) -> Vector3<f64> {
    let v = point - plane_point;
    let d = v.dot(&plane_normal);
    point - d * plane_normal
}

/// Given a point ALREADY PROJECTED INTO THE PLANE, and the three points of a triangle, determine if the point is inside the triangle
/// https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/point_in_triangle.html
pub fn is_point_in_triangle(
    point: Vector3<f64>,
    ta: Vector3<f64>,
    tb: Vector3<f64>,
    tc: Vector3<f64>,
) -> bool {
    // Make point the origin
    let a = ta - point;
    let b = tb - point;
    let c = tc - point;

    // Compute the normals
    let u = b.cross(&c);
    let v = c.cross(&a);
    let w = a.cross(&b);

    if u.dot(&v) < 0.0 {
        return false;
    }
    if u.dot(&w) < 0.0 {
        return false;
    }
    true
}

pub fn vector3_to_array(v: Vector3<f64>) -> [f64; 3] {
    [v.x, v.y, v.z]
}
pub fn quaternion_to_array(q: UnitQuaternion<f64>) -> [f64; 4] {
    [q.i, q.j, q.k, q.w]
}

/// Assuming the three poitns are in COUNTER-CLOCKWISE order, calculate the normal of a (triangular) face
pub fn calculate_face_normal(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> Vector3<f64> {
    let ab = b - a;
    let ac = c - a;
    ab.cross(&ac).normalize()
}

/// Calculate MSE between two matrices
pub fn mse(a: &Matrix3<f64>, b: &Matrix3<f64>) -> f64 {
    let mut sum = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            sum += (a[(i, j)] - b[(i, j)]).powi(2);
        }
    }
    sum / 9.0
}