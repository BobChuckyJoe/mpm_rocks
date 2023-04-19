use crate::{equations::{polar_ru_decomp}, math_utils::mse};
use nalgebra::{Matrix3};

// Granite
pub const GRANITE_DENSITY: f64 = 1463.64 / 4.0; // kg/m^3

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