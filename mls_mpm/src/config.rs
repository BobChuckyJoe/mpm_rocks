pub const SIMULATION_SIZE: f64 = 10.0;
pub const GRID_LENGTH_X: usize = 100;
pub const GRID_LENGTH_Y: usize = 100;
pub const GRID_LENGTH_Z: usize = 100;
pub const GRID_LENGTHS: (usize, usize, usize) = (GRID_LENGTH_X, GRID_LENGTH_Y, GRID_LENGTH_Z);

pub const GRID_SPACING: f64 = SIMULATION_SIZE / GRID_LENGTH as f64;
pub const DELTA_T: f64 = 0.001;
pub const N_PARTICLES: usize = 1000;
pub const N_ITERATIONS: usize = 10000;
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
pub const TIME_TO_SAVE: Option<usize> = Some(N_ITERATIONS - 1);
pub const OUTPUT_GRID_DISTANCES: Option<usize> = TIME_TO_SAVE; // usize is the timestep we want to save
pub const OUTPUT_GRID_DISTANCE_SIGNS: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_VELOCITIES: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_AFFINITIES: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_PARTICLE_DEFORMATION_GRADIENT: Option<usize> = TIME_TO_SAVE;
pub const OUTPUT_GRID_FORCES: Option<usize> = TIME_TO_SAVE;
