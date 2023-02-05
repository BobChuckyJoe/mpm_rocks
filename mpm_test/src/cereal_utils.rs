// Stuff for serializing stuff (JSON). Cereal. Ha. Get it?
use serde_derive::{Deserialize, Serialize};

use crate::types::Particle;

#[derive(Deserialize, Serialize, Debug)]
pub struct Simulation {
    box_size: f64,
    grid_length: usize,
    grid_spacing: f64,
    delta_t: f64,
    num_particles: usize,
    pub num_iterations: usize,
    pub particle_positions: Vec<Vec<[f64; 2]>>,
}

impl Simulation {
    pub fn new (
        box_size: f64,
        grid_length: usize,
        grid_spacing: f64,
        delta_t: f64,
        num_particles: usize,
        num_iterations: usize,
        particle_positions: Vec<Vec<[f64; 2]>>,
    ) -> Simulation {
        Simulation {
            box_size,
            grid_length,
            grid_spacing,
            delta_t,
            num_particles,
            num_iterations,
            particle_positions,
        }
    }

    pub fn add_particle_pos(&mut self, particles: &Vec<Particle>) {
        let mut particle_pos: Vec<[f64; 2]> = Vec::new();
        for p in particles {
            particle_pos.push([p.position.x, p.position.y]);
        }
        self.particle_positions.push(particle_pos);
    }
}