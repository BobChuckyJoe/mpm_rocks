use serde_derive::{Deserialize, Serialize};

use crate::{types::{RigidBody,Particle,Gridcell}, math_utils::{quaternion_to_array, vector3_to_array}};

#[derive(Deserialize, Serialize, Debug)]
pub struct Simulation {
    box_size: f64,
    grid_length: usize,
    grid_spacing: f64,
    delta_t: f64,
    num_particles: usize,
    pub num_iterations: usize,
    pub particle_positions: Vec<Vec<[f64; 3]>>,
    pub rigid_body_positions: Vec<[f64; 3]>,
    pub rigid_body_orientations: Vec<[f64; 4]>,
    pub rigid_body_velocities: Vec<[f64; 3]>,
    pub rigid_body_angular_momentums: Vec<[f64; 3]>,
    pub obj_file_com: [f64; 3], // Center of mass of the obj file
    // Grid things
    pub unsigned_distance_field: Vec<Vec<Vec<Vec<f64>>>>, // At each timestep, the distance of each grid point to the closest rigid body
    pub grid_velocities: Vec<Vec<Vec<Vec<[f64; 3]>>>>,
}

impl Simulation {
    pub fn new(
        box_size: f64,
        grid_length: usize,
        grid_spacing: f64,
        delta_t: f64,
        num_particles: usize,
        num_iterations: usize,
        particle_positions: Vec<Vec<[f64; 3]>>,
        rigid_body_positions: Vec<[f64; 3]>,
        rigid_body_orientations: Vec<[f64; 4]>,
        rigid_body_velocities: Vec<[f64; 3]>,
        rigid_body_angular_momentums: Vec<[f64;3]>,
        obj_file_com: [f64; 3],
        unsigned_distance_field: Vec<Vec<Vec<Vec<f64>>>>,
        grid_velocities: Vec<Vec<Vec<Vec<[f64; 3]>>>>,
    ) -> Simulation {
        Simulation {
            box_size,
            grid_length,
            grid_spacing,
            delta_t,
            num_particles,
            num_iterations,
            particle_positions,
            rigid_body_positions,
            rigid_body_orientations,
            rigid_body_velocities,
            rigid_body_angular_momentums,
            obj_file_com,
            unsigned_distance_field,
            grid_velocities: grid_velocities
        }
    }

    pub fn add_particle_pos(&mut self, particles: &Vec<Particle>) {
        let mut particle_pos: Vec<[f64; 3]> = Vec::new();
        for p in particles {
            particle_pos.push([p.position.x, p.position.y, p.position.z]);
        }
        self.particle_positions.push(particle_pos);
    }
    pub fn add_rigid_body_stuff(&mut self, rb: &RigidBody) {
        self.rigid_body_positions.push(vector3_to_array(rb.position));
        self.rigid_body_orientations.push(quaternion_to_array(rb.orientation));
        self.rigid_body_velocities.push(vector3_to_array(rb.velocity));
        self.rigid_body_angular_momentums.push(vector3_to_array(rb.angular_momentum));
    }

    pub fn add_signed_distance_field(&mut self, grid: &Vec<Vec<Vec<Gridcell>>>) {
        let mut to_ret: Vec<Vec<Vec<f64>>> = Vec::new();
        for i in 0..grid.len() {
            let mut inner: Vec<Vec<f64>> = Vec::new();
            for j in 0..grid[i].len() {
                let mut inner_2: Vec<f64> = Vec::new();
                for k in 0..grid[i][j].len() {
                    inner_2.push(grid[i][j][k].unsigned_distance);
                }
                inner.push(inner_2);
            }
            to_ret.push(inner);
        }
        self.unsigned_distance_field.push(to_ret);
    }

    pub fn add_grid_velocities(&mut self, grid: &Vec<Vec<Vec<Gridcell>>>) {
        let mut to_ret: Vec<Vec<Vec<[f64; 3]>>> = Vec::new();
        for i in 0..grid.len() {
            let mut inner: Vec<Vec<[f64; 3]>> = Vec::new();
            for j in 0..grid[i].len() {
                let mut inner_2: Vec<[f64; 3]> = Vec::new();
                for k in 0..grid[i][j].len() {
                    inner_2.push(vector3_to_array(grid[i][j][k].velocity));
                }
                inner.push(inner_2);
            }
            to_ret.push(inner);
        }
        self.grid_velocities.push(to_ret);
    }
}
