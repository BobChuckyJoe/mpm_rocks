use nalgebra::{Vector2, Matrix2};

#[derive(Copy, Clone)]
pub struct Particle {
    pub position: Vector2<f64>,
    pub velocity: Vector2<f64>,
    pub mass: f64,
    pub density: f64,
    // deformation gradient split into elastic and plastic parts
    pub fe: Matrix2<f64>,
    pub fp: Matrix2<f64>,
}

#[derive(Copy, Clone)]
pub struct Cell {
    pub velocity: Vector2<f64>,
    pub mass: f64,
    pub force: Vector2<f64>, // f_i from equation 6
}

/*
Given a particle, give the corresponding grid cell it is contained in
*/
pub fn to_grid(p: Particle, grid_spacing: f64, grid_length: usize) -> (i64, i64) {
    let x = p.position.x;
    let y = p.position.y;
    let i = (x / grid_spacing).floor() as i64;
    let j = (y / grid_spacing).floor() as i64;
    return (i,j);
}