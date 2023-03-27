use std::sync::Arc;

/// Random helper functions for math things
use nalgebra::{Matrix3, Vector3};

use crate::config::GRID_LENGTH;

// Helper function to iterate through 3x3x3 grid around the starting index
// This is mainly used by rigid particle splatting
pub fn iterate_over_3x3(start: (usize, usize, usize)) -> Vec<(usize, usize, usize)> {
    let mut inds = Vec::new();
    let (i, j, k) = start;
    for di in -1..2 {
        for dj in -1..2 {
            for dk in 0..3 {
                let new_i = i as i32 + di;
                let new_j = j as i32 + dj;
                let new_k = k as i32 + dk;
                if (new_i >= 0 && new_i < GRID_LENGTH as i32)
                    && (new_j >= 0 && new_j < GRID_LENGTH as i32)
                    && (new_k >= 0 && new_k < GRID_LENGTH as i32)
                {
                    inds.push((new_i as usize, new_j as usize, new_k as usize));
                }
                inds.push((start.0 + i, start.1 + j, start.2 + k));
            }
        }
    }
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
