pub const SIMULATION_SIZE: f64 = 5.0;
pub const GRID_LENGTH: usize = 10;
pub const GRID_SPACING: f64 = SIMULATION_SIZE / GRID_LENGTH as f64;
pub const DELTA_T: f64 = 0.001;
pub const N_PARTICLES: usize = 4;
pub const N_ITERATIONS: usize = 10;
pub const BOUNDARY: f64 = 0.05; // Particles this close to the boundary have their velocities zeroed out

//  Values from the 2013 snow paper
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

// Path to the rigid body mesh to load
pub const RIGID_BODY_PATH: &str = "icosahedron.obj";
pub const RIGID_BODY_WEIGHT: f64 = 1.0; // Currently assuming isotropic density
pub const RIGID_BODY_PARTICLES_PER_FACE: usize = 30;

// For penalty force
pub const PENALTY_STIFFNESS: f64 = 1e-5;
pub const BOUNDARY_C: f64 = 0.05; // Used to calculate v tilde (in equation 27 and 28 of MLS MPM paper)

// For output things
pub const OUTPUT_GRID_DISTANCES: Option<usize> = Some(0); // usize is the timestep we want to save
pub const OUTPUT_GRID_DISTANCE_SIGNS: Option<usize> = Some(0);
pub const OUTPUT_GRID_VELOCITIES: Option<usize> = Some(0);
pub const OUTPUT_GRID_AFFINITIES: Option<usize> = Some(0);
