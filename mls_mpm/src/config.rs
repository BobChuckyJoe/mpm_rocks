use nalgebra::Vector3;

use crate::particle_init::uniform_sphere_centered_at_middle;
use crate::types::Particle;

pub const SIMULATION_DIMENSIONS: (f64, f64, f64) = (10.0, 10.0, 10.0);
pub const GRID_LENGTH_X: usize = 30;
pub const GRID_LENGTH_Y: usize = 30;
pub const GRID_LENGTH_Z: usize = 30;
pub const GRID_LENGTHS: (usize, usize, usize) = (GRID_LENGTH_X, GRID_LENGTH_Y, GRID_LENGTH_Z);

pub const GRID_SPACING: f64 = SIMULATION_DIMENSIONS.0 / GRID_LENGTH_X as f64;

pub const DELTA_T: f64 = 0.001;
pub const N_PARTICLES: usize = 1000;
pub const N_ITERATIONS: usize = 50;
pub const BOUNDARY: f64 = 4.0 * GRID_SPACING; // Particles this close to the boundary have their velocities zeroed out
pub const DIMENSIONS: usize = 3;

// Path to the rigid body mesh to load
pub const RIGID_BODY_PATH: &str = "icosahedron.obj";
pub const RIGID_BODY_WEIGHT: f64 = 1.0; // Currently assuming isotropic density
pub const RIGID_BODY_PARTICLES_PER_FACE: usize = 30;

// For penalty force
pub const PENALTY_STIFFNESS: f64 = 1e6;
pub const BOUNDARY_C: f64 = 0.1; // Used to calculate v tilde (in equation 27 and 28 of MLS MPM paper)

// For output things
// pub const TIME_TO_SAVE: Option<usize> = Some(N_ITERATIONS - 1);
pub const TIME_TO_SAVE: Option<usize> = Some(1);
pub const OUTPUT_GRID_DISTANCES: Option<usize> = TIME_TO_SAVE; // usize is the timestep we want to save
pub const OUTPUT_GRID_DISTANCE_SIGNS: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_VELOCITIES: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_AFFINITIES: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_PARTICLE_DEFORMATION_GRADIENT: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_FORCES: Option<usize> = TIME_TO_SAVE;
pub const PRINT_TIMINGS: bool = false;

// Particle initialization
pub const PARTICLE_INIT_FUNC: fn() -> Vec<Particle> = || uniform_sphere_centered_at_middle(1.5, 1000.0);

// Rigid body initialization
pub const RIGID_BODY_INITIAL_POSITION: Vector3<f64> = Vector3::new(5.0, 5.0, 8.0);
pub const RIGID_BODY_INITIAL_VELOCITY: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
pub const RIGID_BODY_INITIAL_ANGULAR_MOMENTUM: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);