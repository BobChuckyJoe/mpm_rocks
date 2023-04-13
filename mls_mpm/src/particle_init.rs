use std::f64::consts::PI;

use nalgebra::{Matrix3, Vector3};
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;

use crate::config::{INITIAL_DENSITY, N_PARTICLES, SIMULATION_SIZE};
use crate::types::Particle;

pub fn uniform_sphere_centered_at_middle(radius: f64, density: f64) -> Vec<Particle> {
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
                // + Vector3::new(
                //     SIMULATION_SIZE / 2.0,
                //     SIMULATION_SIZE / 2.0,
                //     SIMULATION_SIZE / 2.0, // TODO change this back to center
                // ),
                + Vector3::new(5.0, 5.0, 5.0),
            velocity: Vector3::zeros(),
            apic_b: Matrix3::zeros(),
            mass: 1.0,
            density: density,
            deformation_gradient: Matrix3::identity(),
            affinity: false,
            tag: 0,
            particle_distance: 0.0,
            particle_normal: Vector3::zeros(),
        };
        points.push(p);
    }
    for p in points.iter() {
        assert!(p.position.x >= 0.0);
        assert!(p.position.y >= 0.0);
        assert!(p.position.z >= 0.0);
        assert!(p.position.x <= SIMULATION_SIZE);
        assert!(p.position.y <= SIMULATION_SIZE);
        assert!(p.position.z <= SIMULATION_SIZE);
    }
    // Set the mass of each particle according to the density
    let tot_volume = 4.0 / 3.0 * PI * radius.powi(3);
    let tot_mass = tot_volume * density;
    assert!(points.len() != 0);
    let mass_per_particle = tot_mass / points.len() as f64;
    
    for p in points.iter_mut() {
        p.mass = mass_per_particle;
    }
    points
}
