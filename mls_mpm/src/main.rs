mod config;
mod equations;
mod icosahedron;
mod final_sim_config;
mod math_utils;
mod obj_loader;
mod particle_init;
mod serialize;
mod tests;
mod types;
mod material_properties;

use std::collections::{HashSet, HashMap};
use std::sync::Arc;
use std::sync::Mutex;

use nalgebra::{Matrix3, Vector3};
use rayon::prelude::*;
// use tobj::{load_obj};
// use obj::{load_obj, Obj};
use tqdm::tqdm;

use crate::final_sim_config::{
    BOUNDARY, BOUNDARY_C, DELTA_T, GRID_LENGTH_X, GRID_LENGTH_Y, GRID_LENGTH_Z, GRID_SPACING, N_ITERATIONS, N_PARTICLES,
    PENALTY_STIFFNESS, SIMULATION_DIMENSIONS, OUTPUT_GRID_AFFINITIES, OUTPUT_GRID_DISTANCES, OUTPUT_GRID_VELOCITIES, RIGID_BODY_PATH, OUTPUT_GRID_DISTANCE_SIGNS, OUTPUT_PARTICLE_DEFORMATION_GRADIENT, OUTPUT_GRID_FORCES, TIME_TO_SAVE, GRID_LENGTHS
};
use crate::equations::{
    convert_direction_to_world_coords, convert_to_world_coords, convert_world_coords_to_local,
    convert_world_direction_to_local, get_base_grid_ind, grad_weighting_function, velocity_projection,
    grid_cell_ind_to_world_coords, proj_r, weighting_function, calculate_center_of_mass, update_orientation, get_inertia_tensor_world, get_omega,
};
use crate::final_sim_config::is_in_bounds;
use crate::material_properties::{GRANITE_DENSITY, DIRT_DENSITY, partial_psi_partial_f, CRITICAL_COMPRESSION, CRITICAL_STRETCH, neo_hookean_partial_psi_partial_f, project_to_yield_surface, H_0, H_2, H_3, H_1, sand_partial_psi_partial_f, SAND_DENSITY};
use crate::obj_loader::load_rigid_body;
use crate::math_utils::{is_point_in_triangle, iterate_over_3x3, project_point_into_plane, calculate_face_normal};
use crate::serialize::Simulation;
use crate::types::{Gridcell, Particle, RigidBody};

// See the paragraph above eq 177 of siggraph mpm course
// const D_INV: Matrix3<f64> = Matrix3::new(3.0 / GRID_SPACING / GRID_SPACING, 0.0, 0.0,
//                            0.0, 3.0 / GRID_SPACING / GRID_SPACING, 0.0,
//                            0.0, 0.0, 3.0 / GRID_SPACING / GRID_SPACING);
const D_INV: f64 = 3.0 / GRID_SPACING / GRID_SPACING; // Constant since we're using a cubic spline


fn main() {
    let mut sim = Simulation::new(
        SIMULATION_DIMENSIONS,
        GRID_LENGTHS,
        GRID_SPACING,
        DELTA_T,
        N_PARTICLES,
        N_ITERATIONS,
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        [0.0, 0.0, 0.0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    // let mut particles: Vec<Particle> = particle_init::uniform_sphere_centered_at_middle(1.5, SAND_DENSITY);
    let mut particles: Vec<Particle> = final_sim_config::slope_particle_init();
    let mut grid: HashMap<(usize, usize, usize), Mutex<Gridcell>> = HashMap::new();

    println!("Current directory: {:?}", std::env::current_dir());

    let mut rigid_body: RigidBody = load_rigid_body(RIGID_BODY_PATH, GRANITE_DENSITY);
    sim.add_rigid_body_mesh_data(&rigid_body);
    rigid_body.position = final_sim_config::RIGID_BODY_INITIAL_POSITION;
    rigid_body.velocity = final_sim_config::RIGID_BODY_INITIAL_VELOCITY;
    rigid_body.angular_momentum = Vector3::new(0.0, 10000.0, 0.0);
    println!("Orientation: {:?}", rigid_body.orientation);
    println!("Rigid body mass: {}", rigid_body.mass);

    println!("Initialization stuff done!");
    for iteration_num in tqdm(0..N_ITERATIONS) {
        println!("START OF INTERATION {}", iteration_num);
        // println!("Rigid body position: {:?}", rigid_body.position);
        // println!("Rigid body velocity: {:?}", rigid_body.velocity);
        // println!("Rigid body orientation: {:?}", rigid_body.orientation);
        // println!("Rigid body angular momentum: {:?}", rigid_body.angular_momentum);
        // println!("particle 0 velocity: {:?}", particles[0].velocity);
        
        let start = std::time::Instant::now();
        let should_save = iteration_num % 20 == 0;
        // let should_save = true;
        // let should_save = false;
        // if iteration_num % (0.04 / DELTA_T) as usize == 0 {
        match TIME_TO_SAVE {
            Some(x) => {
                if iteration_num == x || should_save {
                    sim.add_particle_pos(&particles);
                    sim.add_rigid_body_stuff(&rigid_body);
                }
            }
            None  => {
                if should_save{
                    // Write the locations every 40 miliseconds, which corresponds to 25 fps
                    sim.add_particle_pos(&particles);
                    sim.add_rigid_body_stuff(&rigid_body);
                }
            }
        }
        println!("Time to save: {:?}", start.elapsed());
        
        // Reset grid
        // for i in 0..GRID_LENGTH {
        //     for j in 0..GRID_LENGTH {
        //         for k in 0..GRID_LENGTH {
        //             grid[i][j][k].reset_values();
        //         }
        //     }
        // }   
        let start = std::time::Instant::now();
        for gridcell in grid.values_mut() {
            gridcell.lock().unwrap().reset_values();
        }
        println!("Time to reset grid: {:?}", start.elapsed());
        
        // Update rigid body velocity and momentum at the very end!
        let tot_change_in_angular_momentum: Arc<Mutex<nalgebra::Matrix<f64, nalgebra::Const<3>, nalgebra::Const<1>, nalgebra::ArrayStorage<f64, 3, 1>>>> = Arc::new(Mutex::new(Vector3::new(0.0, 0.0, 0.0)));
        let tot_change_in_linear_velocity = Arc::new(Mutex::new(Vector3::new(0.0, 0.0, 0.0)));
        
        // Calculate a Hashmap that maps each gridcell to the set of particles that are affected by it
        let start = std::time::Instant::now();
        let grid_to_particles: std::sync::RwLock<HashMap<(usize, usize, usize), Mutex<HashSet<usize>>>> = std::sync::RwLock::new(HashMap::new());
        particles.par_iter().enumerate().for_each(|(ind, p)| {
            let gridcell = get_base_grid_ind(&particles[ind].position, GRID_SPACING);
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let gridcell = (gridcell.0 as isize + dx, gridcell.1 as isize + dy, gridcell.2 as isize + dz);
                        if gridcell.0 < 0 || gridcell.0 >= GRID_LENGTH_X as isize ||
                        gridcell.1 < 0 || gridcell.1 >= GRID_LENGTH_Y as isize ||
                        gridcell.2 < 0 || gridcell.2 >= GRID_LENGTH_Z as isize {
                            continue;
                        }
                        let gridcell = (gridcell.0 as usize, gridcell.1 as usize, gridcell.2 as usize);
                        let gridcell_read = grid_to_particles.read().unwrap();
                        let gridcell_set = gridcell_read.get(&gridcell);
                        match gridcell_set {
                            Some(x) => {
                                x.lock().unwrap().insert(ind);
                            }
                            None => {
                                drop(gridcell_read);
                                let mut new_set = HashSet::new();
                                new_set.insert(ind);
                                grid_to_particles.write().unwrap().insert(gridcell, Mutex::new(new_set));
                            }
                        }
                    }
                }
            }
            let gridcell_read = grid_to_particles.read().unwrap();
            let gridcell_set = gridcell_read.get(&gridcell);
            match gridcell_set {
                Some(x) => {
                    x.lock().unwrap().insert(ind);
                }
                None => {
                    drop(gridcell_read);
                    let mut new_set = HashSet::new();
                    new_set.insert(ind);
                    grid_to_particles.write().unwrap().insert(gridcell, Mutex::new(new_set));
                }
            }
        });
        println!("Time to calculate hashmap: {:?}", start.elapsed());

        // Iterate through the hashmap to make sure every thing exists on the grid
        let start = std::time::Instant::now();
        for key in grid_to_particles.read().unwrap().keys() {
            grid.entry(*key).or_insert(Gridcell::new().into());
        }
        println!("Time to make sure everything exists on the grid: {:?}", start.elapsed());

        let start = std::time::Instant::now();
        // Calculate Colored Distance Field for rigid body
        for ind in 0..rigid_body.rigid_particle_positions.len() {
            // Convert to world_coords
            let p_pos =
                convert_to_world_coords(&rigid_body, rigid_body.rigid_particle_positions[ind]);
            let triangle = rigid_body.faces[rigid_body.rigid_particle_triangles[ind]];
            let p_surface = (
                convert_to_world_coords(&rigid_body, rigid_body.vertices[triangle.0]),
                convert_to_world_coords(&rigid_body, rigid_body.vertices[triangle.1]),
                convert_to_world_coords(&rigid_body, rigid_body.vertices[triangle.2]),
            );
            let rp_normal = calculate_face_normal(p_surface.0, p_surface.1, p_surface.2);
            
            // assert!(rp_normal.dot(&(p_surface.0 - p_pos)) < 1e-9, "{}", rp_normal.dot(&(p_surface.0 - p_pos)));
            // assert!(is_point_in_triangle(p_pos, p_surface.0, p_surface.1, p_surface.2));
                
            // calculate minumum distance
            let start = get_base_grid_ind(&p_pos, GRID_SPACING);
            for (neighbor_i, neighbor_j, neighbor_k) in iterate_over_3x3(start, GRID_LENGTHS) {
                let grid_cell_loc =
                    Vector3::new(neighbor_i as f64, neighbor_j as f64, neighbor_k as f64)
                        * GRID_SPACING;
                // Check if projection of the grid cell onto the plane is valid (if not valid, skip)
                let proj = project_point_into_plane(grid_cell_loc, rp_normal, p_pos);
                // assert!((proj - p_surface.0).dot(&rp_normal).abs() < 1e-6, "dot_prod: {}", (proj - p_surface.0).dot(&rp_normal).abs()); // Check if projection is valid
                if !is_point_in_triangle(proj, p_surface.0, p_surface.1, p_surface.2) {
                    continue;
                }

                let dist = Vector3::new(
                    grid_cell_loc.x - p_pos.x,
                    grid_cell_loc.y - p_pos.y,
                    grid_cell_loc.z - p_pos.z,
                )
                .norm();
                if neighbor_i >= GRID_LENGTH_X || neighbor_j >= GRID_LENGTH_Y || neighbor_k >= GRID_LENGTH_Z {
                    continue;
                }
                let hashmap_ind = (neighbor_i, neighbor_j, neighbor_k);
                if grid.get(&hashmap_ind).is_none() {
                    grid.insert((neighbor_i, neighbor_j, neighbor_k), Gridcell::new().into());
                }
                let mut gridcell = grid.get(&hashmap_ind).unwrap().lock().unwrap();
                if gridcell.unsigned_distance < dist {
                    continue;
                }
                gridcell.affinity = true;
                // If it is the shortest distance, set the distance the be the current one
                gridcell.unsigned_distance = dist;
                // Calculate the sign of the distance (inside or outside)
                // Negative means inside the rigid body
                if (grid_cell_loc - proj).dot(&rp_normal) < 0.0 {
                    gridcell.distance_sign = -1;
                } else {
                    gridcell.distance_sign = 1;
                }
            }
        }
        println!("Time to calculate CDF: {:?}", start.elapsed());
        
        let start = std::time::Instant::now();
        // Calculate particle tags T_{pr}
        // Equation 21 from MLS-MPM paper
        particles.par_iter_mut().for_each(|p| {
                let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
                let mut sum: f64 = 0.0;
                for dx in -2..3 {
                    for dy in -2..3 {
                        for dz in -2..3 {
                            let x: i64 = base_coord.0 as i64 + dx;
                            let y: i64 = base_coord.1 as i64 + dy;
                            let z: i64 = base_coord.2 as i64 + dz;
                            if x < 0
                                || x >= GRID_LENGTH_X as i64
                                || y < 0
                                || y >= GRID_LENGTH_Y as i64
                                || z < 0
                                || z >= GRID_LENGTH_Z as i64
                            {
                                continue;
                            }
                            let x = x as usize;
                            let y = y as usize;
                            let z = z as usize;
                            let hashmap_ind = (x, y, z);
                            let gridcell = grid.get(&hashmap_ind);
                            // If the gridcell doesn't exist, then it doesn't have an affinity
                            if gridcell.is_none() {
                                continue;
                            }
                            let gridcell = gridcell.unwrap().lock().unwrap();
                            // There should be no contribution if the grid cell is too far
                            if !gridcell.affinity {
                                continue;
                            }
                            // TODO According to 5.3.3, there needs to be some persistence in the color of the affinity,
                            // but i'm updating it every timestep. Is this an issue?
                            p.affinity = true;
                            sum += weighting_function(p.position, (x, y, z))
                                * gridcell.unsigned_distance
                                * gridcell.distance_sign as f64;
                        }
                    }
                }
                p.particle_distance = sum;
                p.tag = if sum > 0.0 { 1 } else { -1 };
            }
        );
        println!("Time to calculate particle tags: {:?}", start.elapsed());

        let start = std::time::Instant::now();
        // Calculate particle normals and particle distances from section 5.3.2
        particles.par_iter_mut().for_each(|p| {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            let mut sum: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x: i64 = base_coord.0 as i64 + dx;
                        let y: i64 = base_coord.1 as i64 + dy;
                        let z: i64 = base_coord.2 as i64 + dz;
                        if x < 0
                            || x >= GRID_LENGTH_X as i64
                            || y < 0
                            || y >= GRID_LENGTH_Y as i64
                            || z < 0
                            || z >= GRID_LENGTH_Z as i64
                        {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        let hashmap_ind = (x, y, z);
                        let gridcell = grid.get(&hashmap_ind).unwrap().lock().unwrap();
                        let grad_w = grad_weighting_function(p.position, (x, y, z));
                        // TODO Is this right...? Is there some chain rule thing i'm missing?
                        sum += grad_w * gridcell.unsigned_distance * gridcell.distance_sign as f64;
                        //* (gridcell_pos - p.position);
                    }
                }
            }
            if sum.x == 0.0 && sum.y == 0.0 && sum.z == 0.0 {
                p.particle_normal = Vector3::new(0.0, 0.0, 0.0);
            }
            else {
                p.particle_normal = sum.normalize();
            }
        }
        );
        println!("Time to calculate particle normals: {:?}", start.elapsed());

        // CPIC p2g

        let start = std::time::Instant::now();
        //mass
        // for p in particles.iter() {
        //     let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
        //     for dx in -2..3 {
        //         for dy in -2..3 {
        //             for dz in -2..3 {
        //                 let x: i64 = base_coord.0 as i64 + dx;
        //                 let y: i64 = base_coord.1 as i64 + dy;
        //                 let z: i64 = base_coord.2 as i64 + dz;
        //                 if x < 0
        //                     || x >= GRID_LENGTH_X as i64
        //                     || y < 0
        //                     || y >= GRID_LENGTH_Y as i64
        //                     || z < 0
        //                     || z >= GRID_LENGTH_Z as i64
        //                 {
        //                     continue;
        //                 }
        //                 let x = x as usize;
        //                 let y = y as usize;
        //                 let z = z as usize;
        //                 let hashmap_ind = (x, y, z);
        //                 if grid.get(&hashmap_ind).is_none() {
        //                     grid.insert(hashmap_ind, Gridcell::new());
        //                 }
        //                 // Check compatibility. If they're not nearby, then don't need to check compatibility?
        //                 if grid.get(&hashmap_ind).unwrap().affinity {
        //                     if grid.get(&hashmap_ind).unwrap().distance_sign != p.tag {
        //                         continue;
        //                     }    
        //                 }
        //                 grid.get_mut(&hashmap_ind).unwrap().mass += p.mass * weighting_function(p.position, (x, y, z));
        //             }
        //         }
        //     }
        // }
        grid_to_particles.read().unwrap().par_iter().for_each(|(key, particle_set)| {
            for p in particle_set.lock().unwrap().iter() {
                let particle = particles[*p];
                let mut grid_cell = grid.get(key).unwrap().lock().unwrap();
                grid_cell.mass += particle.mass * weighting_function(particle.position, (key.0, key.1, key.2));                
            }
        }
        );
        println!("Time to splat mass: {:?}", start.elapsed());
        
        let start = std::time::Instant::now();
        //velocity
        // for p in particles.iter() {
        //     let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
        //     for dx in -2..3 {
        //         for dy in -2..3 {
        //             for dz in -2..3 {
        //                 let x: i64 = base_coord.0 as i64 + dx;
        //                 let y: i64 = base_coord.1 as i64 + dy;
        //                 let z: i64 = base_coord.2 as i64 + dz;
        //                 if x < 0
        //                     || x >= GRID_LENGTH_X as i64
        //                     || y < 0
        //                     || y >= GRID_LENGTH_Y as i64
        //                     || z < 0
        //                     || z >= GRID_LENGTH_Z as i64
        //                 {
        //                     continue;
        //                 }
        //                 let x = x as usize;
        //                 let y = y as usize;
        //                 let z = z as usize;
        //                 let hashmap_ind = (x, y, z);
        //                 if grid.get(&hashmap_ind).is_none() {
        //                     grid.insert(hashmap_ind, Gridcell::new().into());
        //                 }
        //                 let mut gridcell = grid.get(&hashmap_ind).unwrap().lock().unwrap();
        //                 let grid_mass = gridcell.mass;
        //                 // Check compatibility
        //                 if gridcell.affinity {
        //                     if gridcell.distance_sign != p.tag {
        //                         // TODO Possibly sus
        //                         // let gridcell_pos = Vector3::new(x as f64, y as f64, z as f64) * GRID_SPACING;
        //                         // gridcell.velocity += p.mass
        //                         //     * weighting_function(p.position, (x, y, z))
        //                         //     * (proj_r(&rigid_body, p.velocity, p.particle_normal, gridcell_pos)
        //                         //     + D_INV * p.apic_b * (gridcell_pos - p.position)) / grid_mass;
        //                         continue;
        //                     }
        //                 }
        //                 if grid_mass == 0.0 {
        //                     continue;
        //                 }
        //                 let gridcell_pos =
        //                     Vector3::new(x as f64, y as f64, z as f64) * GRID_SPACING;
        //                 gridcell.velocity += p.mass
        //                     * weighting_function(p.position, (x, y, z))
        //                     * (p.velocity + D_INV * p.apic_b * (gridcell_pos - p.position)) / grid_mass;
        //             }
        //         }
        //     }
        // }
        grid_to_particles.read().unwrap().par_iter().for_each(|(key, particle_set)| {
            for p in particle_set.lock().unwrap().iter() {
                let particle = particles[*p];
                let mut grid_cell = grid.get(key).unwrap().lock().unwrap();
                let gridcell_pos = Vector3::new(key.0 as f64, key.1 as f64, key.2 as f64) * GRID_SPACING;
                let gridcell_mass = grid_cell.mass;
                if grid_cell.affinity {
                    if grid_cell.distance_sign != particle.tag {
                        // TODO Possibly sus
                        // let gridcell_pos = Vector3::new(x as f64, y as f64, z as f64) * GRID_SPACING;
                        // gridcell.velocity += p.mass
                        //     * weighting_function(p.position, (x, y, z))
                        //     * (proj_r(&rigid_body, p.velocity, p.particle_normal, gridcell_pos)
                        //     + D_INV * p.apic_b * (gridcell_pos - p.position)) / grid_mass;
                        continue;
                    }
                }
                if gridcell_mass == 0.0 {
                    continue;
                }
                grid_cell.velocity += particle.mass * weighting_function(particle.position, (key.0, key.1, key.2))
                * (particle.velocity + D_INV * particle.apic_b * (gridcell_pos - particle.position)) / gridcell_mass;                
            }
        }
        );
        println!("Time to splat velocity: {:?}", start.elapsed());

        // Calculate the initial particle density
        if iteration_num == 0 {
            for p in particles.iter_mut() {
                let mut tot_density = 0.0;
                let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
                for dx in -2..3 {
                    for dy in -2..3 {
                        for dz in -2..3 {
                            let x: i64 = base_coord.0 as i64 + dx;
                            let y: i64 = base_coord.1 as i64 + dy;
                            let z: i64 = base_coord.2 as i64 + dz;
                            if x < 0
                                || x >= GRID_LENGTH_X as i64
                                || y < 0
                                || y >= GRID_LENGTH_Y as i64
                                || z < 0
                                || z >= GRID_LENGTH_Z as i64
                            {
                                continue;
                            }
                            let x = x as usize;
                            let y = y as usize;
                            let z = z as usize;
                            let hashmap_ind = (x, y, z);
                            if grid.get(&hashmap_ind).is_none() {
                                grid.insert(hashmap_ind, Gridcell::new().into());
                            }
                            tot_density += grid.get(&hashmap_ind).unwrap().lock().unwrap().mass
                                * weighting_function(p.position, (x, y, z)) / GRID_SPACING.powi(3);
                        }
                    }
                }
                p.density = tot_density;
            }
        }
        
        let start = std::time::Instant::now();
        // Grid force calculation using eq 18 from MLS-MPM paper
        // for p in particles.iter() {
        //     let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
        //     for dx in -2..3 {
        //         for dy in -2..3 {
        //             for dz in -2..3 {
        //                 let x: i64 = base_coord.0 as i64 + dx;
        //                 let y: i64 = base_coord.1 as i64 + dy;
        //                 let z: i64 = base_coord.2 as i64 + dz;
        //                 if x < 0
        //                     || x >= GRID_LENGTH_X as i64
        //                     || y < 0
        //                     || y >= GRID_LENGTH_Y as i64
        //                     || z < 0
        //                     || z >= GRID_LENGTH_Z as i64
        //                 {
        //                     continue;
        //                 }
        //                 let x = x as usize;
        //                 let y = y as usize;
        //                 let z = z as usize;
        //                 let hashmap_ind = (x, y, z);
        //                 if grid.get(&hashmap_ind).is_none() {
        //                     grid.insert(hashmap_ind, Gridcell::new().into());
        //                 }
                        
        //                 let particle_volume = p.mass / p.density;
                        
        //                 let m_inv = D_INV;
        //                 let partial_psi_partial_f = neo_hookean_partial_psi_partial_f(p.f_e * p.f_p);
        //                 let grid_cell_position =
        //                     Vector3::new(x as f64, y as f64, z as f64) * GRID_SPACING;
        //                 grid.get_mut(&hashmap_ind).unwrap().lock().unwrap().force += -weighting_function(p.position, (x, y, z))
        //                     * particle_volume
        //                     * m_inv
        //                     * partial_psi_partial_f
        //                     * (p.f_e * p.f_p).transpose()
        //                     * (grid_cell_position - p.position);
        //                 // if p.f_e.determinant() < 0.1 {
        //                 //     println!("x y z: {} {} {}", x, y, z);
        //                 //     println!("Big compression bad Det: {}", {p.f_e.determinant()});
        //                 //     println!("Particle position: {:?}", p.position);
        //                 //     println!("Particle deformation gradient: {:?}", p.f_e);
        //                 //     println!("particle volume: {}", particle_volume);
        //                 //     println!("m_inv: {}", m_inv);
        //                 //     println!("partial_psi_partial_f: {}", partial_psi_partial_f);
        //                 //     println!("grid_cell_position: {}", grid_cell_position);
        //                 //     println!("p.f_e * p.f_p: {}", p.f_e * p.f_p);
        //                 //     println!("grid_cell_position - p.position: {}", grid_cell_position - p.position);
        //                 //     let diff = -weighting_function(p.position, (x, y, z))
        //                 //     * particle_volume
        //                 //     * m_inv
        //                 //     * partial_psi_partial_f
        //                 //     * (p.f_e * p.f_p).transpose()
        //                 //     * (grid_cell_position - p.position);
        //                 //     println!("diff: {}", diff);
        //                 //     should_break = true;
        //                 // }
        //             }
        //         }

        //     }
        //     // if (p.f_e.determinant() < 0.1) {
        //     //     panic!();
        //     // }
        // }
        grid_to_particles.read().unwrap().par_iter().for_each(|(key, particle_set)| {
            for p in particle_set.lock().unwrap().iter() {
                let particle = particles[*p];
                let mut grid_cell = grid.get(key).unwrap().lock().unwrap();
                let gridcell_pos = Vector3::new(key.0 as f64, key.1 as f64, key.2 as f64) * GRID_SPACING;
                let gridcell_mass = grid_cell.mass;
                let particle_volume = particle.mass / particle.density;
                        
                let m_inv = D_INV;
                let partial_psi_partial_f = sand_partial_psi_partial_f(particle.f_e * particle.f_p);
                grid_cell.force += -weighting_function(particle.position, (key.0, key.1, key.2))
                    * particle_volume
                    * m_inv
                    * partial_psi_partial_f
                    * (particle.f_e * particle.f_p).transpose()
                    * (gridcell_pos - particle.position);
            }
        }
        );
        println!("Time to update grid forces: {:?}", start.elapsed());
        // println!("Grid force of grid[20][17][9]: {}", grid[20][17][9].force);
        // if should_break {
        //     panic!();
        // }
        // Grid velocity update
        // It says that the grid position should not actually change, but the grid velocity should.
        // See eq 182 of mpm course
        // I'll just do explicit Euler for now
        // for i in 0..GRID_LENGTH {
        //     for j in 0..GRID_LENGTH {
        //         for k in 0..GRID_LENGTH {
        //             if grid[i][j][k].mass == 0.0 {
        //                 continue;
        //             }
        //             let grid_force = grid[i][j][k].force;
        //             let grid_mass = grid[i][j][k].mass;
        //             grid[i][j][k].velocity += DELTA_T * grid_force / grid_mass;
        //             // Gravity
        //             grid[i][j][k].velocity += Vector3::new(0.0, 0.0, -9.8) * DELTA_T;
        //         }
        //     }
        // }
        let start = std::time::Instant::now();
        let mut gridcells_to_visit: HashSet<(usize, usize, usize)> = HashSet::new();
        for p in particles.iter() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x: i64 = base_coord.0 as i64 + dx;
                        let y: i64 = base_coord.1 as i64 + dy;
                        let z: i64 = base_coord.2 as i64 + dz;
                        if x < 0
                            || x >= GRID_LENGTH_X as i64
                            || y < 0
                            || y >= GRID_LENGTH_Y as i64
                            || z < 0
                            || z >= GRID_LENGTH_Z as i64
                        {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        gridcells_to_visit.insert((x, y, z));
                    }
                }
            }
        }
        for (x, y, z) in gridcells_to_visit {
            let hashmap_ind = (x, y, z);
            if grid.get(&hashmap_ind).is_none() {
                grid.insert(hashmap_ind, Gridcell::new().into());
            }
            let mut gridcell = grid.get_mut(&hashmap_ind).unwrap().lock().unwrap();
            if gridcell.mass == 0.0 {
                continue;
            }
            let grid_force = gridcell.force;
            let grid_mass = gridcell.mass;
            gridcell.velocity += DELTA_T * grid_force / grid_mass;
            // Gravity
            gridcell.velocity += Vector3::new(0.0, 0.0, -9.8) * DELTA_T;
        }
        println!("Time to update grid velocities: {:?}", start.elapsed());

        // g2p
        
        let start = std::time::Instant::now();
        particles.par_iter_mut().for_each(|p| {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            // Sum these grid cells in the support of the weighting function
            let mut velocity = Vector3::new(0.0, 0.0, 0.0);
            let mut b_new = Matrix3::<f64>::zeros();
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x = base_coord.0 as i64 + dx;
                        let y = base_coord.1 as i64 + dy;
                        let z = base_coord.2 as i64 + dz;
                        if x < 0
                            || x >= GRID_LENGTH_X as i64
                            || y < 0
                            || y >= GRID_LENGTH_Y as i64
                            || z < 0
                            || z >= GRID_LENGTH_Z as i64
                        {
                            continue;
                        }
            
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        let hashmap_ind = (x, y, z);
                        // All cells should be initialized in the g2p at least
                        // if grid.get(&hashmap_ind).is_none() {
                        //     grid.insert(hashmap_ind, Gridcell::new());
                        // }

                        let gridcell = grid.get(&hashmap_ind).unwrap().lock().unwrap();

                        
                        // Check compatibility
                        if gridcell.distance_sign != p.tag && gridcell.affinity && p.affinity {
                            // Incompatible
                            let v_tilde = proj_r(&rigid_body, p.velocity, p.particle_normal, grid_cell_ind_to_world_coords(x, y, z))
                             + (DELTA_T * BOUNDARY_C * p.particle_normal);
                        

                            velocity += weighting_function(p.position, (x, y, z)) * v_tilde;
                            b_new += weighting_function(p.position, (x, y, z))
                                * v_tilde
                                * (grid_cell_ind_to_world_coords(x, y, z) - p.position).transpose();
                            
                            // TODO This part is sus. Check here for bugs
                            // "Each incompatible grid node applies an impulse"
                            let pr = proj_r(
                                &rigid_body,
                                p.velocity,
                                p.particle_normal,
                                grid_cell_ind_to_world_coords(x, y, z),
                            );
                            let impulse = p.mass * (p.velocity - pr) * weighting_function(p.position, (x,y, z));
                            
                            // Linear portion ezpz
                            *tot_change_in_linear_velocity.lock().unwrap() += impulse / rigid_body.mass;
                            
                            let radius = grid_cell_ind_to_world_coords(x, y, z) - rigid_body.position;
                            *tot_change_in_angular_momentum.lock().unwrap() += get_inertia_tensor_world(&rigid_body).try_inverse().unwrap()
                                * radius.cross(&impulse);
                            // println!("new_angular_momentum {:?}", rigid_body.angular_momentum);
                            // println!("new_omega {:?}", rigid_body.omega);
                        } else {
                            // Compatible
                            velocity +=
                                weighting_function(p.position, (x, y, z)) * gridcell.velocity;
                            b_new += weighting_function(p.position, (x, y, z))
                                * gridcell.velocity
                                * (grid_cell_ind_to_world_coords(x, y, z) - p.position).transpose();
                        }
                    }
                }
            }
            p.velocity = velocity;
            // p.velocity += Vector3::new(0.0, 0.0, -9.8 * DELTA_T);
            p.apic_b = b_new;
            
            // Boundary conditions
            // if p.position.x < BOUNDARY
            //     || p.position.y < BOUNDARY
            //     || p.position.z < BOUNDARY
            //     || p.position.x > SIMULATION_SIZE - BOUNDARY
            //     || p.position.y > SIMULATION_SIZE - BOUNDARY
            //     || p.position.z > SIMULATION_SIZE - BOUNDARY
            // {
            //     p.velocity = Vector3::new(0.0, 0.0, 0.0);
            //     p.apic_b = Matrix3::zeros();
            //     p.position.x = p.position.x.clamp(BOUNDARY, SIMULATION_SIZE - BOUNDARY);
            //     p.position.y = p.position.y.clamp(BOUNDARY, SIMULATION_SIZE - BOUNDARY);
            //     p.position.z = p.position.z.clamp(BOUNDARY, SIMULATION_SIZE - BOUNDARY);
            // }
            if !is_in_bounds(p.position) {
                p.velocity = Vector3::new(0.0, 0.0, 0.0);
                p.apic_b = Matrix3::zeros();
            }
        });
        println!("Time to g2p: {:?}", start.elapsed());
        
        let start: std::time::Instant = std::time::Instant::now();
        // Update particle deformation gradient
        particles.par_iter_mut().for_each(|p| {
            let c_n_plus_1 = p.apic_b * D_INV;
            // First, assume all of the deformation is elastic
            let f_e_hat_new =
                (Matrix3::<f64>::identity() + DELTA_T * c_n_plus_1) * p.f_e;
            let svd = f_e_hat_new.svd(true, true);
            // Sand plasticity
            let (new_singular_vals, case, delta_gamma) = project_to_yield_surface(svd, p.alpha);
            p.f_e = svd.u.unwrap() * new_singular_vals * svd.v_t.unwrap();
            p.f_p = p.f_e.try_inverse().unwrap() * f_e_hat_new * p.f_p;
            // Hardening
            let mut delta_q = 0.0;
            match case {
                1 => {
                    delta_q = 0.0;
                }
                2 => {
                    let epsilon_frob_norm = (new_singular_vals[(0, 0)].ln().powi(2) + new_singular_vals[(1, 1)].ln().powi(2) + new_singular_vals[(2, 2)].ln().powi(2)).sqrt();
                    delta_q = epsilon_frob_norm;
                }
                3 => {
                    delta_q = delta_gamma;
                }
                _ => {
                    panic!("Invalid case");
                }
            }
            p.q = p.q + delta_q;
            let phi_f = H_0 + (H_1 * p.q - H_3) * (-H_2 * p.q).exp();
            p.alpha = (2.0 / 3.0 as f64).sqrt() * (2.0 * phi_f) / (3.0 - phi_f.sin());
            // Push excess deformation into the plastic part
            // let f_e_new_svd = f_e_new.svd(true, true);
            // let mut f_e_singular = Matrix3::<f64>::new(
            //     f_e_new_svd.singular_values[0],
            //     0.0,
            //     0.0,
            //     0.0,
            //     f_e_new_svd.singular_values[1],
            //     0.0,
            //     0.0,
            //     0.0,
            //     f_e_new_svd.singular_values[2],
            // );
            // let u = f_e_new_svd.u.unwrap();
            // let v_t = f_e_new_svd.v_t.unwrap();
            // f_e_singular[(0,0)] = clamp(f_e_singular[(0,0)], 1.0 - CRITICAL_COMPRESSION, 1.0 + CRITICAL_STRETCH);
            // f_e_singular[(1,1)] = clamp(f_e_singular[(0,0)], 1.0 - CRITICAL_COMPRESSION, 1.0 + CRITICAL_STRETCH);
            // f_e_singular[(2,2)] = clamp(f_e_singular[(0,0)], 1.0 - CRITICAL_COMPRESSION, 1.0 + CRITICAL_STRETCH);
            // // Update deformation gradient
            // p.f_e = f_e_new_svd.u.unwrap() * f_e_singular * f_e_new_svd.v_t.unwrap();
            // p.f_p = v_t.transpose() * f_e_singular.try_inverse().unwrap() * u.transpose() * p.f_p;
            // Particle advection
            p.position += DELTA_T * p.velocity;
            }
        );
        println!("Time to update deformation gradients: {:?}", start.elapsed());

        let start = std::time::Instant::now();
        // Penalty impulse
        for p in particles.iter_mut() {
            if p.particle_distance < 0.0 {
                // println!("PENALTY_IMPULE!");
                let penalty_impulse =
                    -PENALTY_STIFFNESS * p.particle_distance * p.particle_normal * DELTA_T;
                // Particle velocity update
                p.velocity += penalty_impulse / p.mass;
                // Rigid body update
                let rb_impulse = -penalty_impulse;
                // Linear portion
                *tot_change_in_linear_velocity.lock().unwrap() += rb_impulse / rigid_body.mass;
                // Angular portion
                let radius = p.position - rigid_body.position;
                *tot_change_in_angular_momentum.lock().unwrap() += rigid_body.inertia_tensor.try_inverse().unwrap()
                    * radius.cross(&rb_impulse);
            }
        }
        println!("Time to penalty impulse: {:?}", start.elapsed());
        
        // Adding gravity here... is there better way to do this?
        *tot_change_in_linear_velocity.lock().unwrap() += Vector3::new(0.0, 0.0, -9.8 * DELTA_T);
        
        // Rigid body velocity/angular momentum update
        rigid_body.velocity += *tot_change_in_linear_velocity.lock().unwrap();
        rigid_body.angular_momentum += *tot_change_in_angular_momentum.lock().unwrap();

        // Boundary condition for rigid body
        if rigid_body.position.x < final_sim_config::BOUNDARY
            || rigid_body.position.y < final_sim_config::BOUNDARY
            || rigid_body.position.z < final_sim_config::BOUNDARY
            || rigid_body.position.x > final_sim_config::SIMULATION_DIMENSIONS.0 - final_sim_config::BOUNDARY
            || rigid_body.position.y > final_sim_config::SIMULATION_DIMENSIONS.1 - final_sim_config::BOUNDARY
            || rigid_body.position.z > final_sim_config::SIMULATION_DIMENSIONS.2 - final_sim_config::BOUNDARY
        {
            rigid_body.velocity = Vector3::new(0.0, 0.0, 0.0);
        }

        // Rigid body advection
        rigid_body.position += DELTA_T * rigid_body.velocity;
        rigid_body.orientation = update_orientation(rigid_body.orientation, get_omega(&rigid_body));

        // Add grid values to sim
        // match OUTPUT_GRID_DISTANCES {
        //     Some(iteration_to_save) => {
        //         if iteration_num == iteration_to_save {
        //             sim.add_signed_distance_field(&grid);
        //         }
        //     }
        //     None => {},
        // }
        // match OUTPUT_GRID_VELOCITIES {
        //     Some(iteration_to_save) => {
        //         if iteration_num == iteration_to_save {
        //             sim.add_grid_velocities(&grid);
        //         }
        //     }
        //     None => {},
        // }
        // match OUTPUT_GRID_AFFINITIES {
        //     Some(iteration_to_save) => {
        //         if iteration_num == iteration_to_save {
        //             sim.add_grid_affinities(&grid);
        //         }
        //     },
        //     None => {},
        // }
        // match OUTPUT_GRID_DISTANCE_SIGNS {
        //     Some(iteration_to_save) => {
        //         if iteration_num == iteration_to_save {
        //             sim.add_grid_distance_signs_from_hashmap(&mut grid, GRID_LENGTHS);
        //         }
        //     },
        //     None => {},
        // }
        // match OUTPUT_PARTICLE_DEFORMATION_GRADIENT {
        //     Some(iteration_to_save) => {
        //         if iteration_num == iteration_to_save {
        //             sim.add_particle_deformation_gradients(&particles);
        //         }
        //     },
        //     None => {},
        // }
        // match OUTPUT_GRID_FORCES {
        //     Some(iteration_to_save) => {
        //         if iteration_num == iteration_to_save {
        //             sim.add_grid_forces(&grid);
        //         }
        //     },
        //     None => {},
        // }
    }
    // Since I'm dropping frames, the total number of "iterations" is different. Need to rename
    sim.num_iterations = sim.particle_positions.len();
    std::fs::write("sim.json", serde_json::to_string_pretty(&sim).unwrap()).unwrap();
}
