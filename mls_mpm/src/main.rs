mod config;
mod equations;
mod icosahedron;
mod math_utils;
mod obj_loader;
mod particle_init;
mod serialize;
mod tests;
mod types;
mod material_properties;

use nalgebra::{Matrix3, Vector3};
// use tobj::{load_obj};
// use obj::{load_obj, Obj};
use tqdm::tqdm;

use crate::config::{
    BOUNDARY, BOUNDARY_C, DELTA_T, GRID_LENGTH, GRID_SPACING, N_ITERATIONS, N_PARTICLES,
    PENALTY_STIFFNESS, SIMULATION_SIZE, OUTPUT_GRID_AFFINITIES, OUTPUT_GRID_DISTANCES, OUTPUT_GRID_VELOCITIES, RIGID_BODY_PATH, OUTPUT_GRID_DISTANCE_SIGNS, OUTPUT_PARTICLE_DEFORMATION_GRADIENT, OUTPUT_GRID_FORCES
};
use crate::equations::{
    convert_direction_to_world_coords, convert_to_world_coords, convert_world_coords_to_local,
    convert_world_direction_to_local, get_base_grid_ind, grad_weighting_function, velocity_projection,
    grid_cell_ind_to_world_coords, partial_psi_partial_f, proj_r, weighting_function, calculate_center_of_mass, update_orientation, get_inertia_tensor_world, get_omega,
};
use crate::icosahedron::create_icosahedron;
use crate::material_properties::{GRANITE_DENSITY, DIRT_DENSITY};
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
        SIMULATION_SIZE,
        GRID_LENGTH,
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
    let mut particles: Vec<Particle> = particle_init::uniform_sphere_centered_at_middle(1.5, DIRT_DENSITY);
    let mut grid: Vec<Vec<Vec<Gridcell>>> =
        vec![vec![vec![Gridcell::new(); GRID_LENGTH]; GRID_LENGTH]; GRID_LENGTH];

    println!("Current directory: {:?}", std::env::current_dir());

    let mut rigid_body: RigidBody = load_rigid_body(RIGID_BODY_PATH, GRANITE_DENSITY);
    sim.add_rigid_body_mesh_data(&rigid_body);
    rigid_body.position = Vector3::new(5.0, 5.0, 8.0);
    rigid_body.angular_momentum = Vector3::new(0.0, 0.0, 0.0);
    println!("Orientation: {:?}", rigid_body.orientation);
    println!("Rigid body mass: {}", rigid_body.mass);

    println!("Initialization stuff done!");
    for iteration_num in tqdm(0..N_ITERATIONS) {
        // println!("START OF INTERATION {}", iteration_num);
        // println!("Rigid body position: {:?}", rigid_body.position);
        // println!("Rigid body velocity: {:?}", rigid_body.velocity);
        // println!("Rigid body orientation: {:?}", rigid_body.orientation);
        // println!("Rigid body angular momentum: {:?}", rigid_body.angular_momentum);
        // println!("particle 0 velocity: {:?}", particles[0].velocity);
        // TODO Remove later
        rigid_body.position = Vector3::new(5.0, 5.0, 8.0);
    
        // if iteration_num % (0.04 / DELTA_T) as usize == 0 {
            // Write the locations every 40 miliseconds, which corresponds to 25 fps
            sim.add_particle_pos(&particles);
            sim.add_rigid_body_stuff(&rigid_body);
        // }
        // Reset grid
        for i in 0..GRID_LENGTH {
            for j in 0..GRID_LENGTH {
                for k in 0..GRID_LENGTH {
                    grid[i][j][k].reset_values();
                }
            }
        }
        // Update rigid body velocity and momentum at the very end!
        let mut tot_change_in_angular_momentum = Vector3::new(0.0, 0.0, 0.0);
        let mut tot_change_in_linear_velocity = Vector3::new(0.0, 0.0, 0.0);
        /*
        The code here is APIC.
        I am doing CPIC for the rigid body stuff, which comes after this
        // APIC p2g mass transfer
        for p in particles.iter() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x = base_coord.0 as i64 + dx;
                        let y = base_coord.1 as i64 + dy;
                        let z = base_coord.2 as i64 + dz;
                        if x < 0 || x >= GRID_LENGTH as i64 || y < 0 || y >= GRID_LENGTH as i64 || z < 0 || z >= GRID_LENGTH as i64 {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        let weight = weighting_function(p.position, (x, y, z));
                        grid[x][y][z].mass += weight * p.mass;
                    }
                }
            }
        }
        // APIC p2g momentum transfer
        for p in particles.iter() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);

            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x = base_coord.0 as i64 + dx;
                        let y = base_coord.1 as i64 + dy;
                        let z = base_coord.2 as i64 + dz;
                        if x < 0 || x >= GRID_LENGTH as i64 || y < 0 || y >= GRID_LENGTH as i64 || z < 0 || z >= GRID_LENGTH as i64 {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        if grid[x][y][z].mass == 0.0 {
                            continue;
                        }

                        let grid_world_coords = Vector3::new(x as f64, y as f64, z as f64) * GRID_SPACING;
                        let weight = weighting_function(p.position, (x, y, z));
                        let mass = grid[x][y][z].mass;
                        grid[x][y][z].velocity += weight * p.mass *
                        (p.velocity + p.apic_b * D_INV * (grid_world_coords - p.position))
                        / mass;
                    }
                }
            }
        }
        // Grid force calculation using eq 18 from MLS-MPM paper
        for p in particles.iter() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x: i64 = base_coord.0 as i64 + dx;
                        let y: i64 = base_coord.1 as i64 + dy;
                        let z: i64 = base_coord.2 as i64 + dz;
                        if x < 0 || x >= GRID_LENGTH as i64 || y < 0 || y >= GRID_LENGTH as i64 || z < 0 || z >= GRID_LENGTH as i64 {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        let particle_volume = p.mass / p.density;
                        let m_inv = D_INV;
                        let partial_psi_partial_f = partial_psi_partial_f(p.deformation_gradient);
                        let grid_cell_position = Vector3::new(x as f64, y as f64, z as f64) * GRID_SPACING;
                        grid[x][y][z].force += - weighting_function(p.position, (x, y, z)) * particle_volume * m_inv * partial_psi_partial_f * p.deformation_gradient.transpose() * (grid_cell_position - p.position);
                    }
                }
            }
        }
        // Grid velocity update
        // TODO I'm not exactly understanding how to implement the semi-implicit Euler integration.
        // The update is done at a previous timestep
        // It says that the grid position should not actually change, but the grid velocity should.
        // See eq 182 of mpm course
        // I'll just do explicit Euler for now
        for i in 0..GRID_LENGTH {
            for j in 0..GRID_LENGTH {
                for k in 0..GRID_LENGTH {
                    if grid[i][j][k].mass == 0.0 {
                        continue;
                    }
                    let grid_force = grid[i][j][k].force;
                    let grid_mass = grid[i][j][k].mass;
                    grid[i][j][k].velocity += DELTA_T * grid_force / grid_mass;
                    grid[i][j][k].velocity += Vector3::new(0.0, 0.0, -9.8) * DELTA_T; // Gravity
                }
            }
        }
        */
        // Calculate Colored Distance Field for rigid body
        for ind in 0..rigid_body.rigid_particle_positions.len() {
            // Convert to world_coords
            let p_pos =
                convert_to_world_coords(&rigid_body, rigid_body.rigid_particle_positions[ind]);
            // println!("rigid_body orientation  {:?}", rigid_body.orientation); TODO
            // println!("{:?}", rigid_body.orientation.to_rotation_matrix());
            // println!("p_pos: {:?}", p_pos);
            // println!("rigid_particle ind {}", ind);
            // println!("p_pos: {:?}", p_pos);
            let triangle = rigid_body.faces[rigid_body.rigid_particle_triangles[ind]];
            let p_surface = (
                convert_to_world_coords(&rigid_body, rigid_body.vertices[triangle.0]),
                convert_to_world_coords(&rigid_body, rigid_body.vertices[triangle.1]),
                convert_to_world_coords(&rigid_body, rigid_body.vertices[triangle.2]),
            );
            // let rp_normal = convert_direction_to_world_coords(
            //     &rigid_body,
            //     rigid_body.rigid_particle_normals[ind],
            // );
            let rp_normal = calculate_face_normal(p_surface.0, p_surface.1, p_surface.2);
            // println!("p_pos: {:?}", p_pos);
            assert!(rp_normal.dot(&(p_surface.0 - p_pos)) < 1e-9, "{}", rp_normal.dot(&(p_surface.0 - p_pos)));
            assert!(is_point_in_triangle(p_pos, p_surface.0, p_surface.1, p_surface.2));
                
            // calculate minumum distance
            let start = get_base_grid_ind(&p_pos, GRID_SPACING);
            for (neighbor_i, neighbor_j, neighbor_k) in iterate_over_3x3(start) {
                let grid_cell_loc =
                    Vector3::new(neighbor_i as f64, neighbor_j as f64, neighbor_k as f64)
                        * GRID_SPACING;
                // Check if projection of the grid cell onto the plane is valid (if not valid, skip)
                let proj = project_point_into_plane(grid_cell_loc, rp_normal, p_pos);
                // assert!(&(p_surface.1 - p_surface.0).cross(&(p_surface.2 - p_surface.0)).dot(&rp_normal).abs() < &1e-6); // Check if the normal is valid
                assert!((proj - p_surface.0).dot(&rp_normal).abs() < 1e-6, "dot_prod: {}", (proj - p_surface.0).dot(&rp_normal).abs()); // Check if projection is valid
                if !is_point_in_triangle(proj, p_surface.0, p_surface.1, p_surface.2) {
                    continue;
                }

                let dist = Vector3::new(
                    grid_cell_loc.x - p_pos.x,
                    grid_cell_loc.y - p_pos.y,
                    grid_cell_loc.z - p_pos.z,
                )
                .norm();
                if neighbor_i >= GRID_LENGTH || neighbor_j >= GRID_LENGTH || neighbor_k >= GRID_LENGTH {
                    continue;
                }
                if grid[neighbor_i][neighbor_j][neighbor_k].unsigned_distance < dist {
                    continue;
                }
                grid[neighbor_i][neighbor_j][neighbor_k].affinity = true;
                // If it is the shortest distance, set the distance the be the current one
                grid[neighbor_i][neighbor_j][neighbor_k].unsigned_distance = dist;
                // Calculate the sign of the distance (inside or outside)
                // Negative means inside the rigid body
                if (grid_cell_loc - proj).dot(&rp_normal) < 0.0 {
                    grid[neighbor_i][neighbor_j][neighbor_k].distance_sign = -1;
                    // println!("grid cell {:?} is now at distance -{}", (neighbor_i, neighbor_j, neighbor_k), dist);
                } else {
                    grid[neighbor_i][neighbor_j][neighbor_k].distance_sign = 1;
                    // println!("grid cell {:?} is now at distance {}", (neighbor_i, neighbor_j, neighbor_k), dist);
                }
            }
        }
        
        // Calculate particle tags T_{pr}
        // Equation 21 from MLS-MPM paper
        for p in particles.iter_mut() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            let mut sum: f64 = 0.0;
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x: i64 = base_coord.0 as i64 + dx;
                        let y: i64 = base_coord.1 as i64 + dy;
                        let z: i64 = base_coord.2 as i64 + dz;
                        if x < 0
                            || x >= GRID_LENGTH as i64
                            || y < 0
                            || y >= GRID_LENGTH as i64
                            || z < 0
                            || z >= GRID_LENGTH as i64
                        {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        let gridcell = grid[x][y][z];
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
        // Calculate particle normals and particle distances from section 5.3.2
        for p in particles.iter_mut() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            let mut sum: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x: i64 = base_coord.0 as i64 + dx;
                        let y: i64 = base_coord.1 as i64 + dy;
                        let z: i64 = base_coord.2 as i64 + dz;
                        if x < 0
                            || x >= GRID_LENGTH as i64
                            || y < 0
                            || y >= GRID_LENGTH as i64
                            || z < 0
                            || z >= GRID_LENGTH as i64
                        {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        let gridcell = grid[x][y][z];
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
        // CPIC p2g

        //mass
        for p in particles.iter() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x: i64 = base_coord.0 as i64 + dx;
                        let y: i64 = base_coord.1 as i64 + dy;
                        let z: i64 = base_coord.2 as i64 + dz;
                        if x < 0
                            || x >= GRID_LENGTH as i64
                            || y < 0
                            || y >= GRID_LENGTH as i64
                            || z < 0
                            || z >= GRID_LENGTH as i64
                        {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        // Check compatibility. If they're not nearby, then don't need to check compatibility?
                        if grid[x][y][z].affinity {
                            if grid[x][y][z].distance_sign != p.tag {
                                continue;
                            }    
                        }
                        grid[x][y][z].mass += p.mass * weighting_function(p.position, (x, y, z));
                    }
                }
            }
        }
        
        //velocity
        for p in particles.iter() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x: i64 = base_coord.0 as i64 + dx;
                        let y: i64 = base_coord.1 as i64 + dy;
                        let z: i64 = base_coord.2 as i64 + dz;
                        if x < 0
                            || x >= GRID_LENGTH as i64
                            || y < 0
                            || y >= GRID_LENGTH as i64
                            || z < 0
                            || z >= GRID_LENGTH as i64
                        {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        let gridcell = &mut grid[x][y][z];
                        let grid_mass = gridcell.mass;
                        // Check compatibility
                        if gridcell.affinity {
                            if gridcell.distance_sign != p.tag {
                                // TODO Possibly sus
                                // let gridcell_pos = Vector3::new(x as f64, y as f64, z as f64) * GRID_SPACING;
                                // gridcell.velocity += p.mass
                                //     * weighting_function(p.position, (x, y, z))
                                //     * (proj_r(&rigid_body, p.velocity, p.particle_normal, gridcell_pos)
                                //     + D_INV * p.apic_b * (gridcell_pos - p.position)) / grid_mass;
                                continue;
                            }
                        }
                        if grid_mass == 0.0 {
                            continue;
                        }
                        let gridcell_pos =
                            Vector3::new(x as f64, y as f64, z as f64) * GRID_SPACING;
                        grid[x][y][z].velocity += p.mass
                            * weighting_function(p.position, (x, y, z))
                            * (p.velocity + D_INV * p.apic_b * (gridcell_pos - p.position)) / grid_mass;
                    }
                }
            }
        }

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
                                || x >= GRID_LENGTH as i64
                                || y < 0
                                || y >= GRID_LENGTH as i64
                                || z < 0
                                || z >= GRID_LENGTH as i64
                            {
                                continue;
                            }
                            let x = x as usize;
                            let y = y as usize;
                            let z = z as usize;
                            tot_density += grid[x][y][z].mass
                                * weighting_function(p.position, (x, y, z)) / GRID_SPACING.powi(3);
                        }
                    }
                }
            }
        }

        // Grid force calculation using eq 18 from MLS-MPM paper
        for p in particles.iter() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x: i64 = base_coord.0 as i64 + dx;
                        let y: i64 = base_coord.1 as i64 + dy;
                        let z: i64 = base_coord.2 as i64 + dz;
                        if x < 0
                            || x >= GRID_LENGTH as i64
                            || y < 0
                            || y >= GRID_LENGTH as i64
                            || z < 0
                            || z >= GRID_LENGTH as i64
                        {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        let particle_volume = p.mass / p.density;
                        let m_inv = D_INV;
                        let partial_psi_partial_f = partial_psi_partial_f(p.deformation_gradient);
                        let grid_cell_position =
                            Vector3::new(x as f64, y as f64, z as f64) * GRID_SPACING;
                        grid[x][y][z].force += -weighting_function(p.position, (x, y, z))
                            * particle_volume
                            * m_inv
                            * partial_psi_partial_f
                            * p.deformation_gradient.transpose()
                            * (grid_cell_position - p.position);
                    }
                }
            }
        }
        // Grid velocity update
        // It says that the grid position should not actually change, but the grid velocity should.
        // See eq 182 of mpm course
        // I'll just do explicit Euler for now
        for i in 0..GRID_LENGTH {
            for j in 0..GRID_LENGTH {
                for k in 0..GRID_LENGTH {
                    if grid[i][j][k].mass == 0.0 {
                        continue;
                    }
                    let grid_force = grid[i][j][k].force;
                    let grid_mass = grid[i][j][k].mass;
                    grid[i][j][k].velocity += DELTA_T * grid_force / grid_mass;
                    // Gravity
                    grid[i][j][k].velocity += Vector3::new(0.0, 0.0, -9.8) * DELTA_T;
                }
            }
        }

        // g2p
        for p in particles.iter_mut() {
            let base_coord = get_base_grid_ind(&p.position, GRID_SPACING);
            // Currently using sticky boundary
            let v_tilde = velocity_projection(&rigid_body, p.position) + (DELTA_T * BOUNDARY_C * p.particle_normal);
            let mut velocity = Vector3::new(0.0, 0.0, 0.0);
            let mut b_new = Matrix3::<f64>::zeros();
            for dx in -2..3 {
                for dy in -2..3 {
                    for dz in -2..3 {
                        let x = base_coord.0 as i64 + dx;
                        let y = base_coord.1 as i64 + dy;
                        let z = base_coord.2 as i64 + dz;
                        if x < 0
                            || x >= GRID_LENGTH as i64
                            || y < 0
                            || y >= GRID_LENGTH as i64
                            || z < 0
                            || z >= GRID_LENGTH as i64
                        {
                            continue;
                        }
                        let x = x as usize;
                        let y = y as usize;
                        let z = z as usize;
                        let gridcell = grid[x][y][z];
                        // Check compatibility
                        if gridcell.distance_sign != p.tag && gridcell.affinity && p.affinity {
                            // Incompatible
                            velocity += weighting_function(p.position, base_coord) * v_tilde;
                            b_new += weighting_function(p.position, base_coord)
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
                            tot_change_in_linear_velocity += impulse / rigid_body.mass;
                            
                            let radius = grid_cell_ind_to_world_coords(x, y, z) - rigid_body.position;
                            tot_change_in_angular_momentum += get_inertia_tensor_world(&rigid_body).try_inverse().unwrap()
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
            if p.position.x < BOUNDARY
                || p.position.y < BOUNDARY
                || p.position.z < BOUNDARY
                || p.position.x > SIMULATION_SIZE - BOUNDARY
                || p.position.y > SIMULATION_SIZE - BOUNDARY
                || p.position.z > SIMULATION_SIZE - BOUNDARY
            {
                p.velocity = Vector3::new(0.0, 0.0, 0.0);
            }
        }
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
                tot_change_in_linear_velocity += rb_impulse / rigid_body.mass;
                // Angular portion
                let radius = p.position - rigid_body.position;
                tot_change_in_angular_momentum += rigid_body.inertia_tensor.try_inverse().unwrap()
                    * radius.cross(&rb_impulse);
            }
        }

        // Update particle deformation gradient
        for p in particles.iter_mut() {
            let c_n_plus_1 = p.apic_b * D_INV;
            let f_new =
                (Matrix3::<f64>::identity() + DELTA_T * c_n_plus_1) * p.deformation_gradient;
            p.deformation_gradient = f_new;
            // TODO Do the clamping thing on the singular values for snow plasticity

            // Particle advection
            p.position += DELTA_T * p.velocity;
        }
        
        // Adding gravity here... is there better way to do this?
        tot_change_in_linear_velocity += Vector3::new(0.0, 0.0, -9.8 * DELTA_T);
        
        // Boundary condition for rigid body
        if rigid_body.position.x < BOUNDARY
            || rigid_body.position.y < BOUNDARY
            || rigid_body.position.z < BOUNDARY
            || rigid_body.position.x > SIMULATION_SIZE - BOUNDARY
            || rigid_body.position.y > SIMULATION_SIZE - BOUNDARY
            || rigid_body.position.z > SIMULATION_SIZE - BOUNDARY
        {
            rigid_body.velocity = Vector3::new(0.0, 0.0, 0.0);
        }
        // Rigid body velocity/angular momentum update
        rigid_body.velocity += tot_change_in_linear_velocity;
        rigid_body.angular_momentum += tot_change_in_angular_momentum;

        // Rigid body advection
        rigid_body.position += DELTA_T * rigid_body.velocity;
        rigid_body.orientation = update_orientation(rigid_body.orientation, get_omega(&rigid_body));

        // Add grid values to sim
        match OUTPUT_GRID_DISTANCES {
            Some(iteration_to_save) => {
                if iteration_num == iteration_to_save {
                    sim.add_signed_distance_field(&grid);
                }
            }
            None => {},
        }
        match OUTPUT_GRID_VELOCITIES {
            Some(iteration_to_save) => {
                if iteration_num == iteration_to_save {
                    sim.add_grid_velocities(&grid);
                }
            }
            None => {},
        }
        match OUTPUT_GRID_AFFINITIES {
            Some(iteration_to_save) => {
                if iteration_num == iteration_to_save {
                    sim.add_grid_affinities(&grid);
                }
            },
            None => {},
        }
        match OUTPUT_GRID_DISTANCE_SIGNS {
            Some(iteration_to_save) => {
                if iteration_num == iteration_to_save {
                    sim.add_grid_distance_signs(&grid);
                }
            },
            None => {},
        }
        match OUTPUT_PARTICLE_DEFORMATION_GRADIENT {
            Some(iteration_to_save) => {
                if iteration_num == iteration_to_save {
                    sim.add_particle_deformation_gradients(&particles);
                }
            },
            None => {},
        }
        match OUTPUT_GRID_FORCES {
            Some(iteration_to_save) => {
                if iteration_num == iteration_to_save {
                    sim.add_grid_forces(&grid);
                }
            },
            None => {},
        }
    }
    // Since I'm dropping frames, the total number of "iterations" is different. Need to rename
    sim.num_iterations = sim.particle_positions.len();
    std::fs::write("sim.json", serde_json::to_string_pretty(&sim).unwrap()).unwrap();
}
