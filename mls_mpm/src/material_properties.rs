use crate::{equations::{polar_ru_decomp}, math_utils::mse, config::DIMENSIONS};
use nalgebra::{Matrix3, SVD, Const};

// Granite
pub const GRANITE_DENSITY: f64 = 1463.64; // kg/m^3

// Dirt
// Values taken from https://www.engineeringtoolbox.com/dirt-mud-densities-d_1727.html
pub const DIRT_DENSITY: f64 = 1220.0;

// Snow
// Values taken from 2013 
pub const CRITICAL_COMPRESSION: f64 = 2.5e-2; // theta_c
pub const CRITICAL_STRETCH: f64 = 7.5e-3; // theta_s
pub const HARDENING_COEFFICIENT: f64 = 10.0; // Epsilon
pub const INITIAL_DENSITY: f64 = 4e2; // rho_0 I think this is the density of the material, not of each material point...?
pub const YOUNGS_MODULUS: f64 = 1.4e5; // E_0
pub const POISSONS_RATIO: f64 = 0.2; // nu
                                     // What are these parameters and what do they mean exactly
pub const MU_0: f64 = YOUNGS_MODULUS / (2.0 * (1.0 + POISSONS_RATIO)); // lame parameter
pub const LAMBDA_0: f64 =
    YOUNGS_MODULUS * POISSONS_RATIO / ((1.0 + POISSONS_RATIO) * (1.0 - 2.0 * POISSONS_RATIO)); // lame parameter

pub fn snow_mew(fp: Matrix3<f64>) -> f64 {
    MU_0 * (HARDENING_COEFFICIENT * ( 1.0 - fp.determinant())).exp()
}

pub fn snow_lambda(fp: Matrix3<f64>) -> f64 {
    LAMBDA_0 * (HARDENING_COEFFICIENT * ( 1.0 - fp.determinant())).exp()
}
pub fn partial_psi_partial_f(f_e: Matrix3<f64>, f_p: Matrix3<f64>) -> Matrix3<f64>{
    let (r, u) = polar_ru_decomp(f_e);
    assert!(mse(&(r * u), &f_e) < 1e-5, "r: {}, u: {}, f_e: {}", r, u, f_e);
    let _res = r.map(|x| assert!(x.is_finite()));
    let _res = u.map(|x| assert!(x.is_finite()));

    // println!("fe: {}", fe);
    // println!("r: {}", r);
    // println!("u: {}", u);
    // let f_e_adjugate = f_e.determinant() * f_e.try_inverse().unwrap();
    // let adj_transpose = f_e_adjugate.transpose();

    let mut partial: Matrix3<f64> = Matrix3::zeros();
    // TODO Check this with eq 1
    // for i in 0..3 {
    //     for j in 0..3 {
    //         partial[(i,j)] = 2.0 * snow_mew(f_p) * (f_e[(i,j)] - r[(i,j)]) 
    //         + snow_lambda(f_p) * (f_e.determinant() - 1.0) * adj_transpose[(i,j)];
    //         let fe_ij = f_e[(i,j)];
    //         assert!(fe_ij.is_finite());
    //         assert!(r[(i,j)].is_finite());
    //         // assert!(mew(fp).is_finite(), "mew: {}, fp: {}", mew(fp), fp);
    //         assert!(snow_lambda(f_p).is_finite());
    //         assert!(f_e.determinant().is_finite());
    //         assert!(partial[(i,j)].is_finite());
    //     }
    // }
    
    // Using the 1st Kirchoff Stress Tensor from the Elasticity notes
    // TODO This isn't exactly right, since the snow_mew and snow_lambda are functions of fp, which the derivative doesn't take into account
    partial = 2.0 * MU_0 * (f_e - r) 
        + LAMBDA_0 * (r.transpose() * (f_e * f_p) - Matrix3::identity()) * r;
    partial
}

pub fn neo_hookean_partial_psi_partial_f(deformation_gradient: Matrix3<f64>) -> Matrix3<f64> {
    // if deformation_gradient.determinant() < 1e-5 {
    //     println!("Determinant is too small");
    //     println!("Determinant: {}", deformation_gradient.determinant());
    // }
    MU_0 * (deformation_gradient - deformation_gradient.transpose().try_inverse().unwrap()) + 
    LAMBDA_0 / 2.0 * deformation_gradient.determinant().powi(2).log(10.0) * deformation_gradient.transpose().try_inverse().unwrap()
}

// Drucker-Prager sand
const SAND_YOUNGS_MODULUS: f64 = 3.5e7;
const SAND_POISSONS_RATIO: f64 = 0.2;
const SAND_MU_0: f64 = SAND_YOUNGS_MODULUS / (2.0 * (1.0 + SAND_POISSONS_RATIO));
const SAND_LAMBDA_0: f64 = SAND_YOUNGS_MODULUS * SAND_POISSONS_RATIO / ((1.0 + SAND_POISSONS_RATIO) * (1.0 - 2.0 * SAND_POISSONS_RATIO)); // lame parameter
pub const H_0: f64 = 35.0;
pub const H_1: f64 = 0.0;
pub const H_2: f64 = 0.2;
pub const H_3: f64 = 10.0;
pub const SAND_DENSITY: f64 = 2200.0;

pub fn sand_partial_psi_partial_f(deformation_gradien: Matrix3<f64>) -> Matrix3<f64>{
    let svd = deformation_gradien.svd(true, true);
    let singular_val_inv = Matrix3::new(
        1.0 / svd.singular_values[0],
        0.0,
        0.0,
        0.0,
        1.0 / svd.singular_values[1],
        0.0,
        0.0,
        0.0,
        1.0 / svd.singular_values[2],
    );
    let ln_signular_val = Matrix3::new(
        svd.singular_values[0].ln(),
        0.0,
        0.0,
        0.0,
        svd.singular_values[1].ln(),
        0.0,
        0.0,
        0.0,
        svd.singular_values[2].ln(),
    );

    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();

    u * (2.0 * SAND_MU_0 * singular_val_inv * ln_signular_val + SAND_LAMBDA_0 * ln_signular_val.trace() * singular_val_inv) * v_t.transpose()
    
}
pub fn project_to_yield_surface(f_e_svd: SVD<f64, Const<3>, Const<3>>, particle_alpha: f64) -> (Matrix3<f64>, usize, f64) {
    let epsilon = Matrix3::new(
        f_e_svd.singular_values[0].ln(), 0.0, 0.0,
        0.0, f_e_svd.singular_values[1].ln(), 0.0,
        0.0, 0.0, f_e_svd.singular_values[2].ln(),
    );
    // Equation 27 from sand paper
    let epsilon_hat = epsilon - epsilon.trace() / DIMENSIONS as f64 * Matrix3::identity();
    let epsilon_hat_frobenius_norm = (epsilon_hat[(0,0)].powi(2) + epsilon_hat[(1,1)].powi(2) + epsilon_hat[(2,2)].powi(2)).sqrt();
    let delta_gamma = epsilon_hat_frobenius_norm +
    (DIMENSIONS as f64 * SAND_LAMBDA_0 + 2.0 * SAND_MU_0) / (2.0 * SAND_MU_0) * epsilon.trace() * particle_alpha;

    if delta_gamma <= 0.0 {
        // Case 1
        (Matrix3::new(f_e_svd.singular_values[0], 0.0, 0.0,
                      0.0, f_e_svd.singular_values[1], 0.0,
                      0.0, 0.0, f_e_svd.singular_values[2]), 1, 0.0)
    }
    else if epsilon_hat_frobenius_norm.abs() < 1e-9 || epsilon.trace() > 0.0 {
        // Case 2
        (Matrix3::identity(), 2, 0.0) 
    }
    else {
        // Case 3
        let mut h_p = epsilon - delta_gamma * epsilon_hat / epsilon_hat_frobenius_norm;
        h_p[(0,0)] = h_p[(0,0)].exp();
        h_p[(1,1)] = h_p[(1,1)].exp();
        h_p[(2,2)] = h_p[(2,2)].exp();
        (h_p, 3, delta_gamma)
    }

}