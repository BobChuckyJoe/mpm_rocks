use rand_chacha::ChaCha8Rng;
use rand::{Rng, SeedableRng};
use nalgebra::{Matrix3, Vector3};

use crate::{types::Particle, material_properties::{SAND_DENSITY, H_0, H_1, H_3, H_2}};

pub const SIMULATION_DIMENSIONS: (f64, f64, f64) = (33.0, 5.0, 33.0);
pub const GRID_SPACING: f64 = 0.1;
pub const GRID_LENGTHS: (usize, usize, usize) = ((SIMULATION_DIMENSIONS.0 as f64 / GRID_SPACING) as usize,
                                                 (SIMULATION_DIMENSIONS.1 as f64 / GRID_SPACING) as usize,
                                                 (SIMULATION_DIMENSIONS.2 as f64 / GRID_SPACING) as usize);
pub const DELTA_T: f64 = 0.0001;
pub const N_PARTICLES: usize = 100000;
pub const N_ITERATIONS: usize = 100000;
pub const SOIL_THICCNESS: f64 = 0.5;
pub const BOUNDARY: f64 = 4.0 * GRID_SPACING; // Particles this close to the boundary have their velocities zeroed out
pub const DIMENSIONS: usize = 3;
pub const BOUNDARY_C: f64 = 0.1; // Used to calculate v tilde (in equation 27 and 28 of MLS MPM paper)
pub const GRID_LENGTH_X: usize = GRID_LENGTHS.0;
pub const GRID_LENGTH_Y: usize = GRID_LENGTHS.1;
pub const GRID_LENGTH_Z: usize = GRID_LENGTHS.2;

pub const PENALTY_STIFFNESS: f64 = 1e6;

// Rigid body initial state
pub const RIGID_BODY_PATH: &str = "icosahedron.obj";
pub const RIGID_BODY_INITIAL_POSITION: Vector3<f64> = Vector3::<f64>::new(3.0, 2.5, 30.0);
pub const RIGID_BODY_INITIAL_VELOCITY: Vector3<f64> = Vector3::<f64>::new(6.0, 0.0, 0.0);
pub const RIGID_BODY_INITIAL_ANGULAR_MOMENTUM: Vector3<f64> = Vector3::<f64>::new(0.0, 0.0, 0.0);

// For output things
pub const TIME_TO_SAVE: Option<usize> = Some(N_ITERATIONS - 1);
pub const OUTPUT_GRID_DISTANCES: Option<usize> = TIME_TO_SAVE; // usize is the timestep we want to save
pub const OUTPUT_GRID_DISTANCE_SIGNS: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_VELOCITIES: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_AFFINITIES: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_PARTICLE_DEFORMATION_GRADIENT: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_FORCES: Option<usize> = TIME_TO_SAVE;
pub const PRINT_TIMINGS: bool = false;

pub const PARTICLE_INIT_FUNC: fn() -> Vec<Particle> = slope_particle_init;

pub const top_origin: Vector3<f64> = Vector3::<f64>::new(14.5, 0.0, 6.0);
pub const top_a_vec: Vector3<f64> = Vector3::<f64>::new(-14.5, 0.0, 24.0);
pub const top_b_vec: Vector3<f64> = Vector3::<f64>::new(0.0, 5.0, 0.0);
// pub const top_plane_normal: Vector3<f64> = top_b_vec.cross(&top_a_vec).normalize();

pub const bot_origin: Vector3<f64> = Vector3::<f64>::new(30.0, 0.0, 0.0);
pub const bot_a_vec: Vector3<f64> = Vector3::<f64>::new(-15.5, 0.0, 6.0);
pub const bot_b_vec: Vector3<f64> = Vector3::<f64>::new(0.0, 5.0, 0.0);
// pub const bot_plane_normal: Vector3<f64> = bot_b_vec.cross(&bot_a_vec).normalize();

pub fn slope_particle_init() -> Vec<Particle>{
    let mut particles_to_ret: Vec<Particle> = vec![];
    let rng = &mut rand_chacha::ChaCha8Rng::seed_from_u64(420); 
    // The upper part of the slope 
    let top_plane_normal: Vector3<f64> = top_b_vec.cross(&top_a_vec).normalize();
    // Total volume
    let volume = top_a_vec.dot(&top_b_vec.cross(&(top_plane_normal * SOIL_THICCNESS))).abs();
    let per_particle_mass = volume * SAND_DENSITY;
    for _ in 0..N_PARTICLES/2 {
        let x = rng.gen_range(0.0..1.0);
        let y = rng.gen_range(0.0..1.0);
        let z = rng.gen_range(0.0..1.0) * SOIL_THICCNESS * top_plane_normal;
        let p = Particle {
            position: top_origin + top_a_vec * x + top_b_vec * y + z + Vector3::new(0.0, 0.0, 1.0),
            velocity: Vector3::zeros(),    
            apic_b: Matrix3::zeros(),
            mass: per_particle_mass,
            density: 0.0, // This is set later
            f_e: Matrix3::identity(),
            f_p: Matrix3::identity(),
            affinity: false,
            tag: 0,
            particle_distance: 0.0,
            particle_normal: Vector3::zeros(),
            q: 0.0,
            alpha: 0.0,
        };
        particles_to_ret.push(p);
        assert!(p.position.x >= 0.0 && p.position.x <= SIMULATION_DIMENSIONS.0, "x: {}", p.position.x);
        assert!(p.position.y >= 0.0 && p.position.y <= SIMULATION_DIMENSIONS.1, "y: {}", p.position.y);
        assert!(p.position.z >= 0.0 && p.position.z <= SIMULATION_DIMENSIONS.2, "z: {}", p.position.z);
    }

    // The lower part of the slope 
    // Total volume
    let bot_plane_normal: Vector3<f64> = bot_b_vec.cross(&bot_a_vec).normalize();
    let volume = bot_a_vec.dot(&bot_b_vec.cross(&(bot_plane_normal * SOIL_THICCNESS))).abs();
    let per_particle_mass = volume * SAND_DENSITY;
    for _ in 0..N_PARTICLES/2 {
        let x = rng.gen_range(0.0..1.0);
        let y = rng.gen_range(0.0..1.0);
        let z = rng.gen_range(0.0..1.0) * SOIL_THICCNESS * bot_plane_normal + Vector3::new(0.0, 0.0, 1.0);
        let p = Particle {
            position: bot_origin + bot_a_vec * x + bot_b_vec * y + z,
            velocity: Vector3::zeros(),    
            apic_b: Matrix3::zeros(),
            mass: per_particle_mass,
            density: 0.0, // This is set later
            f_e: Matrix3::identity(),
            f_p: Matrix3::identity(),
            affinity: false,
            tag: 0,
            particle_distance: 0.0,
            particle_normal: Vector3::zeros(),
            q: 0.0,
            alpha: 0.0,
        };
        particles_to_ret.push(p);

        assert!(p.position.x >= 0.0 && p.position.x <= SIMULATION_DIMENSIONS.0);
        assert!(p.position.y >= 0.0 && p.position.y <= SIMULATION_DIMENSIONS.1);
        assert!(p.position.z >= 0.0 && p.position.z <= SIMULATION_DIMENSIONS.2);
    }
    // Calculate alpha for each particle
    let INIT_Q = 0.0;
    let phi_f = H_0 + (H_1 * INIT_Q - H_3) * (-H_2 * INIT_Q).exp(); 
    let alpha = (2.0 / 3.0 as f64).sqrt() * (2.0 * phi_f.sin()) / (3.0 - phi_f.sin());
    for p in particles_to_ret.iter_mut() {
        p.alpha = alpha;
    }

    particles_to_ret
}

pub fn is_in_bounds(position: Vector3<f64>) -> bool {
    let top_plane_normal: Vector3<f64> = top_b_vec.cross(&top_a_vec).normalize();
    let bot_plane_normal: Vector3<f64> = bot_b_vec.cross(&bot_a_vec).normalize();
    if position.x < 0.0 || position.x > SIMULATION_DIMENSIONS.0 ||
       position.y < 0.0 || position.y > SIMULATION_DIMENSIONS.1 ||
       position.z < 0.0 || position.z > SIMULATION_DIMENSIONS.2 {
        false
    }
    // Check if point is below either plane
    else if (position - top_origin).dot(&top_plane_normal) < 0.0 ||
            (position - bot_origin).dot(&bot_plane_normal) < 0.0 {
        false
    }
    else {
        true
    }
}