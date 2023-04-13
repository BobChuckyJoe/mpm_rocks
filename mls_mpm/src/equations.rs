use nalgebra::{Matrix3, Vector3, UnitQuaternion, Quaternion};
use nalgebra::linalg::SVD;
use rand::Rng;
use rand_chacha::ChaCha8Rng;

use crate::config::{GRID_SPACING, LAMBDA_0, MU_0, DELTA_T};
use crate::types::{BoundaryCondition, RigidBody};

pub fn get_base_grid_ind(p: &Vector3<f64>, grid_spacing: f64) -> (usize, usize, usize) {
    return (
        (p.x / grid_spacing).floor() as usize,
        (p.y / grid_spacing).floor() as usize,
        (p.z / grid_spacing).floor() as usize,
    );
}

pub fn grid_cell_ind_to_world_coords(i: usize, j: usize, k: usize) -> Vector3<f64> {
    let x = i as f64 * GRID_SPACING;
    let y = j as f64 * GRID_SPACING;
    let z = k as f64 * GRID_SPACING;
    Vector3::new(x, y, z)
}

// This will use a slightly different interpolation function
// The 2013 paper uses "dyadic products of one-dimensional cubic B-splines"
// pub fn N(x: f64) ->f64 {
//     if 0.0 <= x.abs() && x.abs() <= 1.0 {
//         0.5 * x.abs().powi(3) - x.abs().powi(2) + 2.0 / 3.0
//     }
//     else if 0.0 <= x.abs() && x.abs() < 2.0 {
//         -1.0 / 6.0 * x.abs().powi(3) + x.abs().powi(2) - 2.0 * x.abs() + 4.0 / 3.0
//     }
//     else {
//         0.0
//     }
// }
// But the siggraph mpm course uses a slightly different one ("but still called cubic B-splines")
pub fn N(x: f64) -> f64 {
    if 0.0 <= x.abs() && x.abs() <= 1.0 {
        0.5 * x.abs().powi(3) - x.abs().powi(2) + 2.0 / 3.0
    } else if 1.0 <= x.abs() && x.abs() < 2.0 {
        1.0 / 6.0 * (2.0 - x.abs()).powi(3)
    } else {
        0.0
    }
}
pub fn N_prime(x: f64) -> f64 {
    if x < 0.0 {
        if 0.0 < x.abs() && x.abs() <= 1.0 {
            -1.5 * x.powi(2) - 2.0 * x
        } else if 1.0 < x.abs() && x.abs() < 2.0 {
            0.5 * (2.0 + x).powi(2)
        } else {
            0.0
        }
    } else {
        if 0.0 < x.abs() && x.abs() <= 1.0 {
            1.5 * x.powi(2) - 2.0 * x
        } else if 1.0 < x.abs() && x.abs() < 2.0 {
            -0.5 * (2.0 - x).powi(2)
        } else {
            0.0
        }
    }
}

pub fn weighting_function(particle_pos: Vector3<f64>, grid_inds: (usize, usize, usize)) -> f64 {
    let x = (particle_pos.x - grid_inds.0 as f64 * GRID_SPACING) / GRID_SPACING;
    let y = (particle_pos.y - grid_inds.1 as f64 * GRID_SPACING) / GRID_SPACING;
    let z = (particle_pos.z - grid_inds.2 as f64 * GRID_SPACING) / GRID_SPACING;
    N(x) * N(y) * N(z)
}
pub fn grad_weighting_function(
    particle_pos: Vector3<f64>,
    grid_inds: (usize, usize, usize),
) -> Vector3<f64> {
    let x = (particle_pos.x - grid_inds.0 as f64 * GRID_SPACING) / GRID_SPACING;
    let y = (particle_pos.y - grid_inds.1 as f64 * GRID_SPACING) / GRID_SPACING;
    let z = (particle_pos.z - grid_inds.2 as f64 * GRID_SPACING) / GRID_SPACING;
    Vector3::new(
        N_prime(x) * N(y) * N(z),
        N(x) * N_prime(y) * N(z),
        N(x) * N(y) * N_prime(z),
    )
}

pub fn polar_ru_decomp(mat: Matrix3<f64>) -> (Matrix3<f64>, Matrix3<f64>) {
    let svd = SVD::try_new(mat, true, true, 1e-12, 0).unwrap();
    let singular_mat = Matrix3::new(
        svd.singular_values[0],
        0.0,
        0.0,
        0.0,
        svd.singular_values[1],
        0.0,
        0.0,
        0.0,
        svd.singular_values[2],
    );
    let diff = svd.recompose().unwrap() - mat;
    // assert!(diff.norm() < 1e-6, "Matrix: {}\n with error {}.\n The diff: {}", mat, diff.norm(), diff);
    // assert!((svd.u.unwrap() * singular_mat * svd.v_t.unwrap() - mat).norm() < 1.0, 
    // "The mat: {:?} is not singular. The error: {}", mat,
    // (svd.u.unwrap() * singular_mat * svd.v_t.unwrap() - mat).norm());
    let r = svd.u.unwrap() * svd.v_t.unwrap();
    let u = svd.v_t.unwrap().transpose() * singular_mat * svd.v_t.unwrap();
    // assert!(
    //     (mat - r * u).norm() < 1.0,
    //     "Polar decomposition failed, actual error: {}",
    //     (mat - r * u).norm()
    // );
    (r, u)
}

// There is a typo in the 88-line mls mpm paper implementation I believe
// This implementation is according to eq 52 of mpm course notes
pub fn partial_psi_partial_f(deformation_gradient: Matrix3<f64>) -> Matrix3<f64> {
    let (r, _u) = polar_ru_decomp(deformation_gradient);
    let j = deformation_gradient.determinant();
    // When to use inverse vs when to use pseudo-inverse?
    
    // Interestingly, this seems much slower (roughly 5x slower...)
    // let deformation_mat_transpose_inv = deformation_gradient.transpose().try_inverse()
    // .expect(format!("Something went wrong {}", deformation_gradient).as_str());

    let deformation_mat_transpose_inv = deformation_gradient.transpose().try_inverse();
    match deformation_mat_transpose_inv {
        Some(inv) => {
            return 2.0 * MU_0 * (deformation_gradient - r) * deformation_gradient.transpose()
                + LAMBDA_0 * (j - 1.0) * j * inv
        }
        None => {
            // Take snapshot of particle gradients
            println!("Deformation inverse failed!");
            println!("The deformation gradient: {}", deformation_gradient);
            panic!();
        }
    }
    // TODO idk if i should use something else
    // let deformation_mat_transpose_inv = deformation_gradient.transpose().pseudo_inverse(1e-12);
    // return 2.0 * MU_0 * (deformation_gradient - r) * deformation_gradient.transpose()
    //             + LAMBDA_0 * (j - 1.0) * j * deformation_mat_transpose_inv.unwrap();
}

pub fn convert_to_world_coords(rb: &RigidBody, particle_pos: Vector3<f64>) -> Vector3<f64> {
    rb.orientation.to_rotation_matrix() * particle_pos + rb.position
}

pub fn convert_world_coords_to_local(rb: &RigidBody, world_pos: Vector3<f64>) -> Vector3<f64> {
    rb.orientation.inverse().to_rotation_matrix() * (world_pos - rb.position)
}

pub fn convert_direction_to_world_coords(
    rb: &RigidBody,
    particle_dir: Vector3<f64>,
) -> Vector3<f64> {
    rb.orientation * particle_dir
}

pub fn convert_world_direction_to_local(rb: &RigidBody, world_dir: Vector3<f64>) -> Vector3<f64> {
    rb.orientation.inverse() * world_dir
}

pub fn get_inertia_tensor_world(rb: &RigidBody) -> Matrix3<f64> {
    rb.orientation.to_rotation_matrix() * rb.inertia_tensor * rb.orientation.to_rotation_matrix().inverse()
}

pub fn get_omega(rb: &RigidBody) -> Vector3<f64> {
    let i_world = get_inertia_tensor_world(rb);
    i_world.try_inverse().unwrap() * rb.angular_momentum
}

/// Calculates the velocity of a rigid body at a given point in world coordinates
pub fn velocity_projection(rb: &RigidBody, point: Vector3<f64>) -> Vector3<f64> {
    rb.velocity + get_omega(rb).cross(&(point - rb.position))
}

/// Used to handle different types of rigid body boundaries
/// See MLS MPM paper equation 25
pub fn proj(v: Vector3<f64>, n: Vector3<f64>, bc: BoundaryCondition, mew: f64) -> Vector3<f64> {
    match bc {
        BoundaryCondition::STICKY => Vector3::zeros(),
        BoundaryCondition::SLIP => v - v.dot(&n) * n,
        BoundaryCondition::SEPARATE => {
            if v.dot(&n) <= 0.0 {
                let v_t = v - v.dot(&n) * n;
                let mut x = v_t.norm() + mew * v.dot(&n);
                if x < 0.0 {
                    x = 0.0;
                }
                x * v_t.normalize()
            } else {
                return v;
            }
        }
    }
}

pub fn proj_r(
    rb: &RigidBody,
    particle_velocity: Vector3<f64>,
    particle_normal: Vector3<f64>,
    grid_cell_world_coords: Vector3<f64>,
) -> Vector3<f64> {
    velocity_projection(rb, grid_cell_world_coords)
        + proj(
            particle_velocity - velocity_projection(rb, grid_cell_world_coords),
            particle_normal,
            BoundaryCondition::STICKY,
            0.05,
        ) // TODO The boundary conditions should be stored in the rigid body. Gonna hard code it for now.
}

/// Calculates the volume of a tetrahedron
/// https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
pub fn calculate_mesh_volume(rb: &RigidBody) -> f64 {
    let mut tot_volume = 0.0;
    for triangle_ind in rb.faces.iter() {
        let p1 = rb.vertices[triangle_ind.0];
        let p2 = rb.vertices[triangle_ind.1];
        let p3 = rb.vertices[triangle_ind.2];
        let v321 = p3.x*p2.y*p1.z;
        let v231 = p2.x*p3.y*p1.z;
        let v312 = p3.x*p1.y*p2.z;
        let v132 = p1.x*p3.y*p2.z;
        let v213 = p2.x*p1.y*p3.z;
        let v123 = p1.x*p2.y*p3.z;
        tot_volume += (1.0/6.0)*(-v321 + v231 + v312 - v132 - v213 + v123);
    }
    assert!(tot_volume > 0.0, "Volume of mesh is negative! ({}", tot_volume);
    tot_volume
}

// Given a mesh, calculate the inertia tensor
// Code taken from https://github.com/erich666/jgt-code/blob/master/Volume_11/Number_2/Kallay2006/Moment_of_Inertia.cpp
// Explaination of the method here: https://stackoverflow.com/questions/809832/how-can-i-compute-the-mass-and-moment-of-inertia-of-a-polyhedron
pub fn calculate_inertia_tensor(rb: &RigidBody) -> Matrix3<f64> {
    let mut inertia_tensor = Matrix3::<f64>::zeros();
    let mut m = 0.0;
    let (mut Cx, mut Cy, mut Cz) = (0.0, 0.0, 0.0); // Centroid
    let (mut xx, mut yy, mut zz, mut yx, mut zx, mut zy) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for triangle_inds in rb.faces.iter() {
        let mut vert_inds = triangle_inds;
        // let vert_a = rb.vertices[triangle_inds.0];
        // let vert_b = rb.vertices[triangle_inds.1];
        // let vert_c = rb.vertices[triangle_inds.2];
        // Check the the triangle is oriented correctly
        // if Matrix3::<f64>::from_columns(&[
        //     vert_a,
        //     vert_b,
        //     vert_c,
        // ]).determinant() < 0.0  {
        //     // Swap two vertices
        //     let temp = vert_inds.1;
        //     vert_inds.1 = vert_inds.2;
        //     vert_inds.2 = temp;
        // }

        let vert_a = rb.vertices[vert_inds.0];
        let vert_b = rb.vertices[vert_inds.1];
        let vert_c = rb.vertices[vert_inds.2];
        // Signed tetrahedron volume
        let tet_vol = vert_a.cross(&vert_b).dot(&vert_c); //TODO: I think the code doesn't divide by 6.0 because it's not needed for the inertia tensor...?

        // Contribution to mass
        m += tet_vol;

        // Contribution to centroid
        let sum = vert_a + vert_b + vert_c;
        Cx += tet_vol * sum.x;
        Cy += tet_vol * sum.y;
        Cz += tet_vol * sum.z;

        // Contribution to inertia tensor
        xx += tet_vol
            * (vert_a.x * vert_a.x + vert_b.x * vert_b.x + vert_c.x * vert_c.x + sum.x * sum.x);
        yy += tet_vol
            * (vert_a.y * vert_a.y + vert_b.y * vert_b.y + vert_c.y * vert_c.y + sum.y * sum.y);
        zz += tet_vol
            * (vert_a.z * vert_a.z + vert_b.z * vert_b.z + vert_c.z * vert_c.z + sum.z * sum.z);
        yx += tet_vol
            * (vert_a.y * vert_a.x + vert_b.y * vert_b.x + vert_c.y * vert_c.x + sum.y * sum.x);
        zx += tet_vol
            * (vert_a.z * vert_a.x + vert_b.z * vert_b.x + vert_c.z * vert_c.x + sum.z * sum.x);
        zy += tet_vol
            * (vert_a.z * vert_a.y + vert_b.z * vert_b.y + vert_c.z * vert_c.y + sum.z * sum.y);
    }

    // GetResults function
    let r = 1.0 / (4.0 * m);
    let (Cx, Cy, Cz) = (Cx * r, Cy * r, Cz * r);
    let m = m / 6.0;

    let r = 1.0 / 120.0;
    let Iyx = yx * r - m * Cy * Cx;
    let Izx = zx * r - m * Cz * Cx;
    let Izy = zy * r - m * Cz * Cy;

    let xx = xx * r - m * Cx * Cx;
    let yy = yy * r - m * Cy * Cy;
    let zz = zz * r - m * Cz * Cz;

    let Ixx = yy + zz;
    let Iyy = zz + xx;
    let Izz = xx + yy;

    inertia_tensor[(0, 0)] = Ixx;
    inertia_tensor[(1, 1)] = Iyy;
    inertia_tensor[(2, 2)] = Izz;
    inertia_tensor[(0, 1)] = Iyx;
    inertia_tensor[(1, 0)] = Iyx;
    inertia_tensor[(0, 2)] = Izx;
    inertia_tensor[(2, 0)] = Izx;
    inertia_tensor[(1, 2)] = Izy;
    inertia_tensor[(2, 1)] = Izy;

    inertia_tensor
}

pub fn calculate_center_of_mass(rb: &RigidBody) -> Vector3<f64> {
    let mut center_of_mass = Vector3::<f64>::zeros();
    let mut m = 0.0;
    for triangle_inds in rb.faces.iter() {
        let vert_inds = triangle_inds;
        // let vert_a = rb.vertices[triangle_inds.0];
        // let vert_b = rb.vertices[triangle_inds.1];
        // let vert_c = rb.vertices[triangle_inds.2];
        // Check the the triangle is oriented correctly
        // if Matrix3::<f64>::from_columns(&[
        //     vert_a,
        //     vert_b,
        //     vert_c,
        // ]).determinant() < 0.0  {
        //     // Swap two vertices
        //     let temp = vert_inds.1;
        //     vert_inds.1 = vert_inds.2;
        //     vert_inds.2 = temp;
        // }

        let vert_a = rb.vertices[vert_inds.0];
        let vert_b = rb.vertices[vert_inds.1];
        let vert_c = rb.vertices[vert_inds.2];
        // Signed tetrahedron volume
        let tet_vol = vert_a.cross(&vert_b).dot(&vert_c); //TODO: I think the code doesn't divide by 6.0 because it's not needed for the inertia tensor...?

        // Contribution to mass
        m += tet_vol;

        // Contribution to centroid
        let sum = vert_a + vert_b + vert_c;
        center_of_mass += tet_vol * sum;
    }
    center_of_mass / (4.0 * m)
}

// Using the "reflection" method
// https://blogs.sas.com/content/iml/2020/10/19/random-points-in-triangle.html
pub fn generate_random_point_on_triangle(
    p1: Vector3<f64>,
    p2: Vector3<f64>,
    p3: Vector3<f64>,
    rng: &mut ChaCha8Rng, // For repeatability
) -> Vector3<f64> {
    let a = p2 - p1;
    let b = p3 - p1;
    let mut u1 = rng.gen_range(0.0..1.0);
    let mut u2 = rng.gen_range(0.0..1.0);
    if u1 + u2 > 1.0 {
        u1 = 1.0 - u1;
        u2 = 1.0 - u2;
    }
    let w = u1 * a + u2 * b;
    w + p1
}


pub fn update_orientation(old_q: UnitQuaternion<f64>, omega: Vector3<f64>) -> UnitQuaternion<f64> {
    let omega = omega * DELTA_T;
    if omega.norm() == 0.0 {
        return old_q;
    }
    let rest = omega / omega.norm() * (omega.norm() / 2.0).sin();
    let q = UnitQuaternion::from_quaternion(Quaternion::new(
        (0.5 * omega.norm()).cos(),
        rest.x,
        rest.y,
        rest.z,
    ));
    q * old_q
}