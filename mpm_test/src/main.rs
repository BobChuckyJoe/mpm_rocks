use nalgebra::{Vector2, Matrix2};
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use tqdm::tqdm;

use crate::{types::to_grid, equations::{sigma, grad_weighting_func, CRITICAL_COMPRESSION, CRITICAL_STRETCH}};

mod cereal_utils;
mod equations;
mod types;


const BOX_SIZE: f64 = 10.0; // The size of the box in METERS
const GRID_LENGTH: usize = 1000; // The size of the grid in # of cells
const GRID_SPACING: f64 = BOX_SIZE / GRID_LENGTH as f64; // The spacing between cells in meters
const DELTA_T: f64 = 0.001 ; // Time step in seconds
const N_PARTICLES: usize = 128;
const N_ITERATIONS: usize = 20000;

// // Constituency model parameters from paper
// // TODO This is in 2D so maybe the dimensions aren't right...
// const CRITICAL_COMPRESSION: f64 = 2.5e-2; // theta_c
// const CRITICAL_STRETCH: f64 = 7.5e-3; // theta_s
// const HARDENING_COEFFICIENT: f64 = 10.0; // Epsilon
// const INITIAL_DENSITY: f64 = 4e2; // rho_0
// const YOUNGS_MODULUS: f64 = 1.4e5; // E_0
// const POISSONS_RATIO: f64 = 0.2; // nu


fn main() {
    // Object to write to disk
    let mut sim = cereal_utils::Simulation::new(
        BOX_SIZE,
        GRID_LENGTH,
        GRID_SPACING,
        DELTA_T,
        N_PARTICLES,
        N_ITERATIONS,
        vec![],
    );
    // Initialize simulation space and grid
    let mut GRID: Vec<Vec<types::Cell>> = vec!
        [vec![types::Cell { velocity: Vector2::new(0.0, 0.0), mass: 0.0, force: Vector2::new(0.0, 0.0) }; GRID_LENGTH]; GRID_LENGTH];
    

    // Initialize particles
    // TODO put this on the heap
    let mut PARTICLES: Vec<types::Particle> = vec![types::Particle { 
        position: Vector2::new(0.0, 0.0),
        velocity: Vector2::new(0.0, 0.0),
        mass: 1.0,
        density: 0.0,
        fe: Matrix2::identity(),
        fp: Matrix2::identity(), 
    }
        ; N_PARTICLES];

    // Particle to Grid function things
    let N = |x: f64| {
        if 0.0 <= x.abs() && x.abs() <= 1.0 {
            0.5 * x.abs().powi(3) - x.abs().powi(2) + 2.0 / 3.0
        }
        else if 0.0 <= x.abs() && x.abs() < 2.0 {
            -1.0 / 6.0 * x.abs().powi(3) + x.abs().powi(2) - 2.0 * x.abs() + 4.0 / 3.0
        }
        else {
            0.0
        }
    };
    // My main interpretation of these functions is that they are used to interpolate the particles onto the grid
    let part_to_grid = | grid_ind_i: usize, grid_ind_j: usize, x_particle: f64, y_particle: f64 | {
       return N(1.0/GRID_SPACING * (x_particle - grid_ind_i as f64 * GRID_SPACING)) *
              N(1.0/GRID_SPACING * (y_particle - grid_ind_j as f64 * GRID_SPACING));
    };

    let mut rng = ChaCha8Rng::seed_from_u64(1038);
    // Set initial positions TODO: Make this better or change into different things
    for p in 0..N_PARTICLES {
        // Use rng to set initial positions
        let dx = rng.gen_range(-1.0..1.0);
        let dy = rng.gen_range(-1.0..1.0);
        PARTICLES[p].position.x = BOX_SIZE / 2.0 + dx;
        PARTICLES[p].position.y = BOX_SIZE / 2.0 + dy;

        // Check that particles are inside the grid
        assert!(PARTICLES[p].position.x < BOX_SIZE && PARTICLES[p].position.x >= 0.0);
        assert!(PARTICLES[p].position.y < BOX_SIZE && PARTICLES[p].position.y >= 0.0);
    }


    println!("Initialization stuff done!");
    // exit(0);

    // Main iteration loop
    for curr_timestep in tqdm(0..N_ITERATIONS) {
        // Write down calculations
        // Only write down particle position such that every frame is 1/60th of a second
        if (curr_timestep % (0.016 / DELTA_T).round() as usize == 0) {
            sim.add_particle_pos(&PARTICLES);
        }
        
        // Clear grid
        for x in 0..GRID_LENGTH {
            for y in 0..GRID_LENGTH {
                GRID[x][y].velocity.x = 0.0;
                GRID[x][y].velocity.y = 0.0;
                GRID[x][y].mass = 0.0;
                GRID[x][y].force.x = 0.0;
                GRID[x][y].force.y = 0.0;
            }
        }

        // 1. "Rasterize particle data to the grid"
        // Transfer mass to grid
        for p in PARTICLES.iter() {
            let inds = to_grid(*p, GRID_SPACING, GRID_LENGTH);
            for dx in -1..2  {
                for dy in -1..2 {
                    let new_ind = (inds.0 as i32 + dx, inds.1 as i32 + dy);
                    if new_ind.0 < 0 || new_ind.1 < 0 || new_ind.0 >= GRID_LENGTH as i32 || new_ind.1 >= GRID_LENGTH as i32 {
                        continue;
                    }
                    let weight: f64 = part_to_grid(new_ind.0 as usize, new_ind.1 as usize, p.position.x, p.position.y);
                    GRID[new_ind.0 as usize][new_ind.1 as usize].mass += p.mass * weight;
                }
            }
        }
        // Transfer velocity to grid
        for p in PARTICLES.iter() {
            let inds = to_grid(*p, GRID_SPACING, GRID_LENGTH);
            for dx in -1..2  {
                for dy in -1..2 {
                    let new_ind = (inds.0 as i32 + dx, inds.1 as i32 + dy);
                    if new_ind.0 < 0 || new_ind.1 < 0 || new_ind.0 >= GRID_LENGTH as i32 || new_ind.1 >= GRID_LENGTH as i32 {
                        continue;
                    }
                    let weight: f64 = part_to_grid(new_ind.0 as usize, new_ind.1 as usize, p.position.x, p.position.y);
                    let change_in_velocity = p.velocity * weight * p.mass / GRID[new_ind.0 as usize][new_ind.1 as usize].mass;
                    assert!(GRID[new_ind.0 as usize][new_ind.1 as usize].mass != 0.0);
                    GRID[new_ind.0 as usize][new_ind.1 as usize].velocity += change_in_velocity;
                }
            }
        }
        
        // 2. "Compute particle volume and densities" FIRST TIMESTEP ONLY
        // What is this step for exactly...?
        if curr_timestep == 0 {
            for p in PARTICLES.iter_mut() {
                for i in 0..GRID_LENGTH {
                    for j in 0..GRID_LENGTH {
                        p.density += GRID[i][j].mass * part_to_grid(i, j, p.position.x, p.position.y) / GRID_SPACING.powi(3);
                    }
                }
            }
        }
        
        // 3. "Compute grid forces" with eq 6
        for i in 0..PARTICLES.len() {
            let p = PARTICLES[i];
            let inds = to_grid(p, GRID_SPACING, GRID_LENGTH);
            for dx in -1..2  {
                for dy in -1..2 {
                    let new_ind = (inds.0 as i32 + dx, inds.1 as i32 + dy);
                    if new_ind.0 < 0 || new_ind.1 < 0 || new_ind.0 >= GRID_LENGTH as i32 || new_ind.1 >= GRID_LENGTH as i32 {
                        continue;
                    }
                    // println!("Looking at particle {}", i);
                    GRID[new_ind.0 as usize][new_ind.1 as usize].force -= p.mass / p.density * sigma(p) *
                    grad_weighting_func(new_ind.0 as usize, new_ind.1 as usize, p.position.x, p.position.y, GRID_SPACING); 
                    
                    
                    // Debugging checks
                    let _res = grad_weighting_func(new_ind.0 as usize, new_ind.1 as usize, p.position.x, p.position.y, GRID_SPACING)
                    .map(|x| assert!(x.is_finite())); 
                    let _res = sigma(p).map(|x| assert!(x.is_finite()));
                    assert!(p.density.is_finite());
                    assert!(p.density != 0.0);
                    assert!(p.mass.is_finite());
                    assert!(GRID[new_ind.0 as usize][new_ind.1 as usize].force.x.is_finite());
                    assert!(GRID[new_ind.0 as usize][new_ind.1 as usize].force.y.is_finite());
                }
            }
        }
        
        // 4. "Update velocities on grid to v_i*"
        for x in 0..GRID_LENGTH {
            for y in 0..GRID_LENGTH {
                let inverse_mass = 1.0;
                if (GRID[x][y].mass == 0.0) {
                    
                }
                else {
                    let change_in_velocity = GRID[x][y].force * (1.0 / GRID[x][y].mass) * DELTA_T;
                    let _res = change_in_velocity.map(|x| assert!(x.is_finite())); 
                    GRID[x][y].velocity += change_in_velocity;
                }
                // Add gravity
                GRID[x][y].velocity.x -= 9.8 * DELTA_T;
                // Boundary conditions. TODO: Why is it set to 0 and not -v? 
                if x < 2 || x > GRID_LENGTH - 3 {
                    GRID[x][y].velocity.x = 0.0;
                }
                if y < 2 || y > GRID_LENGTH - 3 {
                    GRID[x][y].velocity.y = 0.0;
                }
            }
        }

        // 5. "Calculate for particle collisions"

        // 6. "Solve the linear system"
        // TODO right now i'm just gonna use the explicit time integration and let v^{n+1} = v^*
        
        // 7. "Update deformation gradient"
        // See list number 7
        for p in PARTICLES.iter_mut() {
            let inds = to_grid(*p, GRID_SPACING, GRID_LENGTH);
            let mut grad_v_p = Matrix2::<f64>::zeros();
            for dx in -1..2 {
                for dy in -1..2 {
                    let new_inds = (inds.0 as i64 + dx, inds.1 as i64 + dy);
                    if new_inds.0 < 0 || new_inds.1 < 0 || new_inds.0 >= GRID_LENGTH as i64 || new_inds.1 >= GRID_LENGTH as i64 {
                        continue;
                    }
                    grad_v_p += GRID[new_inds.0 as usize][new_inds.1 as usize].velocity * grad_weighting_func(new_inds.0 as usize, new_inds.1 as usize,
                         p.position.x, p.position.y, GRID_SPACING).transpose();
                    let _res = GRID[new_inds.0 as usize][new_inds.1 as usize].velocity.map(|x| assert!(x.is_finite()));
                    let _res = p.velocity.map(|x| assert!(x.is_finite()));
                }
            }
            // These updates are defined in section 7
            let fe_hat = (Matrix2::<f64>::identity() + DELTA_T * grad_v_p) * p.fe; // TODO grad_v_p sometimes sus
            let _res = fe_hat.map(|x| assert!(f64::is_finite(x), "{}", grad_v_p));
            // Equation 11
            let f_n_plus_1 = fe_hat * p.fp;

            // "Push" the part of F_e that exceeds the critical deformation threshold to F_p
            let fe_svd = fe_hat.svd(true, true);
            let u = fe_svd.u.unwrap();
            let v_t = fe_svd.v_t.unwrap();
            // Clamp the singular values [1-theta_c, 1+theta_s]
            let clamped_svd = fe_svd.singular_values.map(|x| {
                if x < 1.0 - CRITICAL_COMPRESSION {
                    1.0 - CRITICAL_COMPRESSION
                } else if x > 1.0 + CRITICAL_STRETCH {
                    1.0 + CRITICAL_STRETCH
                } else {
                    x
                }
            });
            // There is prob a better function that can convert this into a diagonal matrix
            let clamped_svd = Matrix2::<f64>::new(*clamped_svd.get(0).unwrap(), 0.0,
                                                                                      0.0, *clamped_svd.get(1).unwrap());
            // Update deformation gradients (eq 12)
            p.fe = u * clamped_svd * v_t;
            let _res = p.fe.map(|x| assert!(f64::is_finite(x), "fe is nan. clamped_svd: {:?}, u: {:?}, v_t: {:?}", clamped_svd, u, v_t));
            // Is the inversion of clamped_svd going to be a problem? Could there be 1 singular value only?
            p.fp = v_t.transpose() * clamped_svd.try_inverse().unwrap() * u.transpose() * f_n_plus_1;
            let _res = p.fp.map(|x| assert!(f64::is_finite(x)));
            assert!((f_n_plus_1 - p.fe * p.fp).norm() < 1e-6);
        }

        // 8. Update particle velocities 
        const alpha: f64 = 0.95;
        for p in PARTICLES.iter_mut() {
            let inds = to_grid(*p, GRID_SPACING, GRID_LENGTH);
            let mut v_pic = Vector2::<f64>::zeros();
            let mut v_flip = p.velocity.clone();
            for dx in -1..2 {
                for dy in -1..2 {
                    let new_inds = (inds.0 as i64 + dx, inds.1 as i64 + dy);
                    if new_inds.0 < 0 || new_inds.1 < 0 || new_inds.0 >= GRID_LENGTH as i64 || new_inds.1 >= GRID_LENGTH as i64 {
                        continue;
                    }
                    v_pic += GRID[new_inds.0 as usize][new_inds.1 as usize].velocity * part_to_grid(new_inds.0 as usize, new_inds.1 as usize, p.position.x, p.position.y);
                    v_flip += (GRID[new_inds.0 as usize][new_inds.1 as usize].velocity - p.velocity) * part_to_grid(new_inds.0 as usize, new_inds.1 as usize, p.position.x, p.position.y);
                }
            }
            // Update velocity
            p.velocity = (1.0 - alpha) * v_pic + alpha * v_flip;
        }

        // 9. "Particle-based body collisions"
        // TODO

        // 10. "Update particle positions"
        for p in PARTICLES.iter_mut() {
            p.position += p.velocity * DELTA_T;
            // Clamp the position to within the range of the grid
            p.position.x = p.position.x.max(0.0).min(GRID_LENGTH as f64 * GRID_SPACING);
            p.position.y = p.position.y.max(0.0).min(GRID_LENGTH as f64 * GRID_SPACING);
        }
        // Keep track of particle deformation gradient
        // println!("Particle 50 fe: {:?}", PARTICLES[50].fe);
        // println!("Particle 50 fp: {:?}", PARTICLES[50].fp);

    }

    std::fs::write("sim.json", serde_json::to_string_pretty(&sim).unwrap()).unwrap();
}
