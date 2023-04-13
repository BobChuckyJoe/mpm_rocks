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
    pub particle_deformation_gradient_x: Vec<Vec<[f64; 3]>>,
    pub particle_deformation_gradient_y: Vec<Vec<[f64; 3]>>,
    pub particle_deformation_gradient_z: Vec<Vec<[f64; 3]>>,
    pub rigid_body_positions: Vec<[f64; 3]>,
    pub rigid_body_orientations: Vec<[f64; 4]>,
    pub rigid_body_velocities: Vec<[f64; 3]>,
    pub rigid_body_angular_momentums: Vec<[f64; 3]>,
    // Rigid body geometry things, doesn't change over time
    pub rigid_body_vertices: Vec<[f64; 3]>,
    pub rigid_body_triangles: Vec<[usize; 3]>,
    pub obj_file_com: [f64; 3], // Center of mass of the obj file
    pub rigid_particle_positions: Vec<[f64; 3]>,
    pub rigid_particle_triangles: Vec<usize>,
    // Grid things
    pub unsigned_distance_field: Vec<Vec<Vec<Vec<f64>>>>, // At each timestep, the distance of each grid point to the closest rigid body
    pub grid_distance_signs: Vec<Vec<Vec<Vec<i32>>>>,
    pub grid_velocities: Vec<Vec<Vec<Vec<[f64; 3]>>>>,
    pub grid_affinities: Vec<Vec<Vec<Vec<bool>>>>,
    pub grid_forces: Vec<Vec<Vec<Vec<[f64; 3]>>>>,
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
        particle_deformation_gradient_x: Vec<Vec<[f64; 3]>>,
        particle_deformation_gradient_y: Vec<Vec<[f64; 3]>>,
        particle_deformation_gradient_z: Vec<Vec<[f64; 3]>>,
        rigid_body_positions: Vec<[f64; 3]>,
        rigid_body_orientations: Vec<[f64; 4]>,
        rigid_body_velocities: Vec<[f64; 3]>,
        rigid_body_angular_momentums: Vec<[f64;3]>,
        rigid_body_vertices: Vec<[f64; 3]>,
        rigid_body_triangles: Vec<[usize; 3]>,
        obj_file_com: [f64; 3],
        rigid_particle_positions: Vec<[f64; 3]>,
        rigid_particle_triangles: Vec<usize>,
        unsigned_distance_field: Vec<Vec<Vec<Vec<f64>>>>,
        grid_distance_signs: Vec<Vec<Vec<Vec<i32>>>>,
        grid_velocities: Vec<Vec<Vec<Vec<[f64; 3]>>>>,
        grid_affinities: Vec<Vec<Vec<Vec<bool>>>>,
        grid_forces: Vec<Vec<Vec<Vec<[f64; 3]>>>>,
    ) -> Simulation {
        Simulation {
            box_size,
            grid_length,
            grid_spacing,
            delta_t,
            num_particles,
            num_iterations,
            particle_positions,
            particle_deformation_gradient_x,
            particle_deformation_gradient_y,
            particle_deformation_gradient_z,
            rigid_body_positions,
            rigid_body_orientations,
            rigid_body_velocities,
            rigid_body_angular_momentums,
            rigid_body_vertices,
            rigid_body_triangles,
            obj_file_com,
            rigid_particle_positions,
            rigid_particle_triangles,
            unsigned_distance_field,
            grid_distance_signs,
            grid_velocities,
            grid_affinities,
            grid_forces,
        }
    }

    pub fn add_particle_pos(&mut self, particles: &Vec<Particle>) {
        let mut particle_pos: Vec<[f64; 3]> = Vec::new();
        for p in particles {
            particle_pos.push([p.position.x, p.position.y, p.position.z]);
        }
        self.particle_positions.push(particle_pos);
    }

    pub fn add_particle_deformation_gradients(&mut self, particles: &Vec<Particle>) {
        let mut def_x: Vec<[f64; 3]> = Vec::new();
        let mut def_y = Vec::new();
        let mut def_z = Vec::new();
        for p in particles {
            let deformation_gradient = p.f_e * p.f_p;
            let x = [deformation_gradient[(0,0)], deformation_gradient[(1,0)], deformation_gradient[(2,0)]];
            def_x.push(x);
            def_y.push([
                deformation_gradient[(0,1)],
                deformation_gradient[(1,1)],
                deformation_gradient[(2,1)],
            ]);
            def_z.push([
                deformation_gradient[(0,2)],
                deformation_gradient[(1,2)],
                deformation_gradient[(2,2)],
            ]);
        }
        self.particle_deformation_gradient_x.push(def_x);
        self.particle_deformation_gradient_y.push(def_y);
        self.particle_deformation_gradient_z.push(def_z);
    }

    pub fn add_rigid_body_stuff(&mut self, rb: &RigidBody) {
        self.rigid_body_positions.push(vector3_to_array(rb.position));
        self.rigid_body_orientations.push(quaternion_to_array(rb.orientation));
        self.rigid_body_velocities.push(vector3_to_array(rb.velocity));
        self.rigid_body_angular_momentums.push(vector3_to_array(rb.angular_momentum));
    }
    pub fn add_rigid_body_mesh_data(&mut self, rb: &RigidBody) {
        self.rigid_body_vertices = rb.vertices.iter().map(|v| vector3_to_array(*v)).collect();
        self.rigid_body_triangles = rb.faces.iter().map(|f| [f.0, f.1, f.2]).collect();
        self.obj_file_com = vector3_to_array(rb.obj_file_com);
        self.rigid_particle_positions = rb.rigid_particle_positions.iter().map(|v| vector3_to_array(*v)).collect();
        self.rigid_particle_triangles = rb.rigid_particle_triangles.clone();
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
    pub fn add_grid_affinities(&mut self, grid: &Vec<Vec<Vec<Gridcell>>>) {
        let mut to_ret: Vec<Vec<Vec<bool>>> = Vec::new();
        for i in 0..grid.len() {
            let mut inner: Vec<Vec<bool>> = Vec::new();
            for j in 0..grid[i].len() {
                let mut inner_2: Vec<bool> = Vec::new();
                for k in 0..grid[i][j].len() {
                    inner_2.push(grid[i][j][k].affinity);
                }
                inner.push(inner_2);
            }
            to_ret.push(inner);
        }
        self.grid_affinities.push(to_ret);
    }
    pub fn add_grid_distance_signs(&mut self, grid: &Vec<Vec<Vec<Gridcell>>>) {
        let mut to_ret: Vec<Vec<Vec<i32>>> = Vec::new();
        for i in 0..grid.len() {
            let mut inner: Vec<Vec<i32>> = Vec::new();
            for j in 0..grid[i].len() {
                let mut inner_2: Vec<i32> = Vec::new();
                for k in 0..grid[i][j].len() {
                    inner_2.push(grid[i][j][k].distance_sign);
                }
                inner.push(inner_2);
            }
            to_ret.push(inner);
        }
        self.grid_distance_signs.push(to_ret);
    }

    pub fn add_grid_forces(&mut self, grid: &Vec<Vec<Vec<Gridcell>>>) {
        let mut to_ret: Vec<Vec<Vec<[f64; 3]>>> = Vec::new();
        for i in 0..grid.len() {
            let mut inner: Vec<Vec<[f64; 3]>> = Vec::new();
            for j in 0..grid[i].len() {
                let mut inner_2: Vec<[f64; 3]> = Vec::new();
                for k in 0..grid[i][j].len() {
                    inner_2.push(vector3_to_array(grid[i][j][k].force));
                }
                inner.push(inner_2);
            }
            to_ret.push(inner);
        }
        self.grid_forces.push(to_ret);
    }
}
