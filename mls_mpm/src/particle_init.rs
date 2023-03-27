use nalgebra::{Matrix3, Vector3};
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;

use crate::config::{INITIAL_DENSITY, N_PARTICLES, SIMULATION_SIZE};
use crate::types::Particle;

pub fn uniform_sphere_centered_at_middle(radius: f64) -> Vec<Particle> {
    const SEED: u64 = 420;
    let mut rng = ChaCha8Rng::seed_from_u64(SEED);
    // Sample from unit ball using acceptance-rejection sampling
    let mut points: Vec<Particle> = vec![];
    while points.len() < N_PARTICLES {
        let x = rng.gen_range(-1.0..1.0);
        let y = rng.gen_range(-1.0..1.0);
        let z = rng.gen_range(-1.0..1.0);
        if x * x + y * y + z * z > 1.0 {
            continue;
        }

        let p = Particle {
            position: Vector3::new(x, y, z) * radius
                + Vector3::new(
                    SIMULATION_SIZE / 2.0,
                    SIMULATION_SIZE / 2.0,
                    SIMULATION_SIZE / 2.0,
                ),
            velocity: Vector3::zeros(),
            apic_b: Matrix3::zeros(),
            mass: 1.0,
            density: INITIAL_DENSITY, //TODO I don't think this is correct...
            deformation_gradient: Matrix3::identity(),
            tag: 0,
            particle_distance: 0.0,
            particle_normal: Vector3::zeros(),
        };
        points.push(p);
    }
    points
}
