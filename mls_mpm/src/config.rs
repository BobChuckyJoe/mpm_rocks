pub const SIMULATION_SIZE: f64 = 10.0;
pub const GRID_LENGTH: usize = 30;
pub const GRID_SPACING: f64 = SIMULATION_SIZE / GRID_LENGTH as f64;
pub const DELTA_T: f64 = 0.001;
pub const N_PARTICLES: usize = 500;
pub const N_ITERATIONS: usize = 2000;
pub const BOUNDARY: f64 = 0.5; // Particles this close to the boundary have their velocities zeroed out

// Path to the rigid body mesh to load
pub const RIGID_BODY_PATH: &str = "icosahedron.obj";
pub const RIGID_BODY_WEIGHT: f64 = 1.0; // Currently assuming isotropic density
pub const RIGID_BODY_PARTICLES_PER_FACE: usize = 30;

// For penalty force
pub const PENALTY_STIFFNESS: f64 = 1e6;
pub const BOUNDARY_C: f64 = 0.1; // Used to calculate v tilde (in equation 27 and 28 of MLS MPM paper)

// For output things
pub const TIME_TO_SAVE: usize = 1140;
pub const OUTPUT_GRID_DISTANCES: Option<usize> = Some(TIME_TO_SAVE); // usize is the timestep we want to save
pub const OUTPUT_GRID_DISTANCE_SIGNS: Option<usize> = Some(TIME_TO_SAVE);
pub const OUTPUT_GRID_VELOCITIES: Option<usize> = Some(TIME_TO_SAVE);
pub const OUTPUT_GRID_AFFINITIES: Option<usize> = Some(TIME_TO_SAVE);
pub const OUTPUT_PARTICLE_DEFORMATION_GRADIENT: Option<usize> = Some(TIME_TO_SAVE);
pub const OUTPUT_GRID_FORCES: Option<usize> = Some(TIME_TO_SAVE);
