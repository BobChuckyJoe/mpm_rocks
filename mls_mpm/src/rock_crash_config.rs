use rand_chacha::ChaCha8Rng;
use rand::{Rng, SeedableRng};
use nalgebra::{Matrix3, Vector3};

use crate::{types::Particle, material_properties::{SAND_DENSITY, H_0, H_1, H_3, H_2}};

pub const SIMULATION_DIMENSIONS: (f64, f64, f64) = (10.0, 10.0, 10.0);
pub const GRID_SPACING: f64 = 0.1;
pub const GRID_LENGTHS: (usize, usize, usize) = ((SIMULATION_DIMENSIONS.0 as f64 / GRID_SPACING) as usize,
                                                 (SIMULATION_DIMENSIONS.1 as f64 / GRID_SPACING) as usize,
                                                 (SIMULATION_DIMENSIONS.2 as f64 / GRID_SPACING) as usize);
pub const DELTA_T: f64 = 0.001;
pub const N_PARTICLES: usize = 1000;
pub const N_ITERATIONS: usize = 100;
pub const SOIL_THICCNESS: f64 = 3.0;
pub const BOUNDARY: f64 = 4.0 * GRID_SPACING; // Particles this close to the boundary have their velocities zeroed out
pub const DIMENSIONS: usize = 3;
pub const BOUNDARY_C: f64 = 0.1; // Used to calculate v tilde (in equation 27 and 28 of MLS MPM paper)
pub const GRID_LENGTH_X: usize = GRID_LENGTHS.0;
pub const GRID_LENGTH_Y: usize = GRID_LENGTHS.1;
pub const GRID_LENGTH_Z: usize = GRID_LENGTHS.2;

pub const PENALTY_STIFFNESS: f64 = 1e6;

// Rigid body initial state
pub const RIGID_BODY_PATH: &str = "icosahedron.obj";
pub const RIGID_BODY_INITIAL_POSITION: Vector3<f64> = Vector3::<f64>::new(5.0, 5.0, 8.0);
pub const RIGID_BODY_INITIAL_VELOCITY: Vector3<f64> = Vector3::<f64>::new(0.0, 0.0, 0.0);
pub const RIGID_BODY_INITIAL_ANGULAR_MOMENTUM: Vector3<f64> = Vector3::<f64>::new(0.0, 100.0, 0.0);

// For output things
pub const TIME_TO_SAVE: Option<usize> = Some(0);
pub const OUTPUT_GRID_DISTANCES: Option<usize> = TIME_TO_SAVE; // usize is the timestep we want to save
pub const OUTPUT_GRID_DISTANCE_SIGNS: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_VELOCITIES: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_AFFINITIES: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_PARTICLE_DEFORMATION_GRADIENT: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_FORCES: Option<usize> = TIME_TO_SAVE;
pub const PRINT_TIMINGS: bool = false;
pub const FILE_OUTPUT_DIR: &str = "rock_crash_small";
// Give time for sand to settle
pub fn should_wait(iteration: usize) -> bool {
    iteration < 100000
}

pub const PARTICLE_INIT_FUNC: fn() -> Vec<Particle> = ground_particle_init;

/// Create a 3m deep thing of sand, and then slam rock into it
pub fn ground_particle_init() -> Vec<Particle>{
    let mut particles_to_ret: Vec<Particle> = vec![];
    let rng = &mut rand_chacha::ChaCha8Rng::seed_from_u64(420); 
    // Total volume
    let volume = 5.0 * 5.0 * 3.0;
    let per_particle_mass = volume * SAND_DENSITY;
    for _ in 0..N_PARTICLES {
        let x = rng.gen_range(2.5..7.5);
        let y = rng.gen_range(2.5..7.5);
        let z = rng.gen_range(0.0..1.0) * SOIL_THICCNESS;
        let p = Particle {
            position: Vector3::new(x, y, z),
            velocity: Vector3::zeros(),    
            apic_b: Matrix3::zeros(),
            mass: per_particle_mass,
            density: SAND_DENSITY,
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

    // Calculate alpha for each particle
    let INIT_Q = 0.0;
    let phi_f = H_0 + (H_1 * INIT_Q - H_3) * (-H_2 * INIT_Q).exp(); 
    let alpha = (2.0 / 3.0 as f64).sqrt() * (2.0 * phi_f.sin()) / (3.0 - phi_f.sin());
    for p in particles_to_ret.iter_mut() {
        p.alpha = alpha;
    }

    particles_to_ret
}
