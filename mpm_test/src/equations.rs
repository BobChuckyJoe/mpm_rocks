use nalgebra::{Vector2, Matrix2};
use crate::types::{Particle};

// Constituency model parameters from paper
// TODO This is in 2D so maybe the dimensions aren't right...
pub const CRITICAL_COMPRESSION: f64 = 2.5e-2; // theta_c
pub const CRITICAL_STRETCH: f64 = 7.5e-3; // theta_s
const HARDENING_COEFFICIENT: f64 = 10.0; // Epsilon
const INITIAL_DENSITY: f64 = 4e2; // rho_0 I think this is the density of the material, not of each material point...?
const YOUNGS_MODULUS: f64 = 1.4e5; // E_0
const POISSONS_RATIO: f64 = 0.2; // nu
// What are these parameters and what do they mean exactly
const mew_nought: f64 = YOUNGS_MODULUS / (2.0 * (1.0 + POISSONS_RATIO)); // lame parameter
const lambda_nought: f64 = YOUNGS_MODULUS * POISSONS_RATIO / ((1.0 + POISSONS_RATIO) * (1.0 - 2.0 * POISSONS_RATIO)); // lame parameter
// Lame (lah-may?) parameters
pub fn mew(fp: Matrix2<f64>) -> f64 {
    mew_nought * (HARDENING_COEFFICIENT * ( 1.0 - fp.determinant())).exp()
}

pub fn lambda(fp: Matrix2<f64>) -> f64 {
    lambda_nought * (HARDENING_COEFFICIENT * ( 1.0 - fp.determinant())).exp()
}

// Implementation of Psi(F_E, F_P)
pub fn energy_density_function(fe: Matrix2<f64>, fp: Matrix2<f64>) -> f64 {
    // TODO I'm using QR decomp right here which might be different from polar decomp
    let polar_decomp = fe.polar();
    println!("{}", polar_decomp.0 * polar_decomp.0.transpose());
    println!("{}", polar_decomp.1 * polar_decomp.1.transpose());
    
    let QR_decomp = fe.qr();
    mew(fp) * frobenius_norm(fe - QR_decomp.q()).powi(2) + 0.5 * lambda(fp) * (fe.determinant() - 1.0).powi(2)
}

/*
    Returns R and U
*/
pub fn polar_ru_decomp(mat: Matrix2<f64>) -> (Matrix2<f64>, Matrix2<f64>) {
    let u_squared = mat.transpose() * mat;
    let mut diag = u_squared.symmetric_eigen();
    
    let diff = (diag.recompose() - u_squared).norm();
    assert!( diff < 1e-6, "Mat: {}, Diff: {}", mat, diff);

    diag.eigenvalues = diag.eigenvalues.map(|x| x.max(0.0).sqrt());
    let u = diag.recompose();
    let r = mat * u.try_inverse().unwrap();
    (r, u)
}

// TODO Not sure if this is how you do this partial (partial psi / partial F_E)
pub fn partial_psi_partial_fe(fe: Matrix2<f64>, fp: Matrix2<f64>) -> Matrix2<f64> {
    let (r, u) = polar_ru_decomp(fe);
    assert!((r * u - fe).norm() < 1e-5);
    let _res = r.map(|x| assert!(x.is_finite()));
    let _res = u.map(|x| assert!(x.is_finite()));

    let mut partial: Matrix2<f64> = Matrix2::zeros();
    // TODO Check this with eq 1
    for i in 0..2 {
        for j in 0..2 {
            partial[(i,j)] = 2.0 * mew(fp) * (fe[(i,j)] - r[(i,j)]) 
            + lambda(fp) * (fe.determinant() - 1.0);
            let fe_ij = fe[(i,j)];
            assert!(fe_ij.is_finite());
            assert!(r[(i,j)].is_finite());
            assert!(mew(fp).is_finite(), "mew: {}, fp: {}", mew(fp), fp);
            assert!(lambda(fp).is_finite());
            assert!(fe.determinant().is_finite());
            assert!(partial[(i,j)].is_finite());
        }
    }

    partial
}

// The weighting function on page 4
pub fn n(x: f64) -> f64 {
    if 0.0 <= x.abs() && x.abs() < 1.0 {
        0.5 * x.abs().powi(3) - x.powi(2) + 2.0/3.0
    }
    else if 1.0 <= x.abs() && x.abs() < 2.0 {
        -1.0/6.0 * x.abs().powi(3) + x.powi(2) - 2.0 * x.abs() + 4.0/3.0
    }
    else {
        0.0
    }
}

pub fn n_prime(x: f64) -> f64 {
    if 0.0 <= x.abs() && x.abs() < 1.0 {
        1.5 * x.powi(2) - 2.0 * x
    } else if 1.0 <= x.abs() && x.abs() < 2.0 {
        -0.5 * x.powi(2) + 2.0 * x - 2.0
    } else {
        0.0
    }
}

pub fn grad_weighting_func(grid_ind_x: usize, grid_ind_y: usize, p_x: f64, p_y: f64, grid_spacing: f64) -> Vector2<f64> {
    let x = n_prime(1.0 / grid_spacing * (p_x - grid_ind_x as f64 * grid_spacing)) * 1.0 / grid_spacing *
     n(1.0 / grid_spacing * (p_y - grid_ind_y as f64 * grid_spacing)); 
    let y = n_prime(1.0 / grid_spacing * (p_y - grid_ind_y as f64 * grid_spacing)) * 1.0 / grid_spacing *
     n(1.0 / grid_spacing * (p_x - grid_ind_x as f64 * grid_spacing));
    
    assert!(x.is_finite());
    assert!(y.is_finite());
    Vector2::new(x, y)
}

// Cauchy stress (eq 6)
pub fn sigma(p: Particle) -> Matrix2<f64> {
    let deformation_matrix = p.fe * p.fp;
    assert!(deformation_matrix.determinant().is_finite());
    let _res = p.fe.map(|x| assert!(x.is_finite()));
    let _res = p.fp.map(|x| assert!(x.is_finite()));
    let _res = partial_psi_partial_fe(p.fe, p.fp).map(|x| assert!(x.is_finite()));
    assert_non_nan(p.fe);
    1.0 / deformation_matrix.determinant() * partial_psi_partial_fe(p.fe, p.fp) * p.fe.transpose()
    
}

pub fn frobenius_norm(matrix: Matrix2<f64>) -> f64 {
    let mut total = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            total += matrix.get((i,j)).unwrap().powi(2);
        }
    }

    total.sqrt()
}

pub fn assert_non_nan(matrix: Matrix2<f64>) {
    for i in 0..2 {
        for j in 0..2 {
            assert!(matrix.get((i,j)).unwrap().is_finite());
        }
    }
}