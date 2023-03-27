use nalgebra::{Matrix3, UnitQuaternion, Vector3};
use tobj::{self, Mesh};

use crate::equations::calculate_center_of_mass;
use crate::types::RigidBody;

pub fn load_rigid_body(path: &str) -> RigidBody {
    let (models, _) = tobj::load_obj(path, &tobj::LoadOptions::default()).unwrap();
    let model = models[0].clone();
    let mesh = model.mesh;

    let position = Vector3::<f64>::zeros();
    let velocity = Vector3::<f64>::zeros();
    let orientation = UnitQuaternion::<f64>::identity();
    let omega = Vector3::<f64>::zeros();
    let mass = 1.0; // TODO: Calculate mass
    let vertices = calculate_vertices(&mesh);
    let faces = calculate_faces(&mesh);
    let vertex_normals = calculate_vertex_normals(&mesh);
    assert!(vertices.len() == vertex_normals.len());
    
    let mut rb = RigidBody {
        position,
        velocity,
        orientation,
        omega,
        mass,
        vertices,
        faces,
        vertex_normals,
        moment_of_inertia: Matrix3::<f64>::identity(), // Added later
        rigid_particle_positions: Vec::new(), // Added later
        rigid_particle_triangles: Vec::new(), // Added later
        rigid_particle_normals: Vec::new(), // Added later
    };
    // Calculate center of mass, and shift eveything if needed
    let com = calculate_center_of_mass(&rb);
    print!("Center of mass: {:?}", com);
    for v in &mut rb.vertices {
        *v -= com;
    }
    // Create rigid particles

    rb
    
}

fn calculate_vertices(mesh: &Mesh) -> Vec<Vector3<f64>> {
    let mut vertices = Vec::new();
    for ind in 0..mesh.positions.len() / 3 {
        let v = Vector3::<f64>::new(
            mesh.positions[3 * ind] as f64,
            mesh.positions[3 * ind + 1] as f64,
            mesh.positions[3 * ind + 2] as f64,
        );
        vertices.push(v);
    }
    vertices
}

fn calculate_faces(mesh: &Mesh) -> Vec<(usize, usize, usize)> {
    let mut faces = Vec::new();
    for ind in 0..mesh.indices.len() / 3 {
        let f = (
            mesh.indices[3 * ind] as usize,
            mesh.indices[3 * ind + 1] as usize,
            mesh.indices[3 * ind + 2] as usize,
        );
        faces.push(f);
    }
    faces
}

fn calculate_vertex_normals(mesh: &Mesh) -> Vec<Vector3<f64>> {
    let mut vertex_normals = Vec::new();
    for ind in 0..mesh.normals.len() / 3 {
        let v = Vector3::<f64>::new(
            mesh.normals[3 * ind] as f64,
            mesh.normals[3 * ind + 1] as f64,
            mesh.normals[3 * ind + 2] as f64,
        );
        vertex_normals.push(v);
    }
    vertex_normals
}

fn create_rigid_particles(rb: &mut RigidBody, num_particles_per_face: usize) {
    for face in rb.faces.iter() {
        let v1 = rb.vertices[face.0];
        let v2 = rb.vertices[face.1];
        let v3 = rb.vertices[face.2];
        
    }
}