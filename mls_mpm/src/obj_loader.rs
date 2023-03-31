use nalgebra::{Matrix3, UnitQuaternion, Vector3};
use tobj::{self, Mesh};

use crate::config::RIGID_BODY_PARTICLES_PER_FACE;
use crate::equations::{calculate_inertia_tensor, calculate_center_of_mass, generate_random_point_on_triangle};
use crate::types::RigidBody;

pub fn load_rigid_body(path: &str) -> RigidBody {
    let (models, _) = tobj::load_obj(path, &tobj::LoadOptions::default()).unwrap();
    let model = models[0].clone();
    let mesh = model.mesh;

    let position = Vector3::<f64>::zeros();
    let velocity = Vector3::<f64>::zeros();
    let orientation = UnitQuaternion::<f64>::identity();
    let angular_momentum = Vector3::<f64>::zeros();
    let mass = 1.0; // TODO: Calculate mass
    let vertices = calculate_vertices(&mesh);
    let faces = calculate_faces(&mesh);
    let vertex_normals = calculate_vertex_normals(&mesh);
    // TODO There are more vertex normals (vn) than vertices... weird
    // assert!(vertices.len() == vertex_normals.len(), "vertices.len() = {}, vertex_normals.len() = {}", vertices.len(), vertex_normals.len());
    
    let mut rb = RigidBody {
        position,
        velocity,
        orientation,
        angular_momentum,
        mass,
        vertices,
        faces,
        vertex_normals,
        inertia_tensor: Matrix3::<f64>::identity(), // Added later
        rigid_particle_positions: Vec::new(), // Added later
        rigid_particle_triangles: Vec::new(), // Added later
        rigid_particle_normals: Vec::new(), // TODO Don't use this, this is wrong... calculate face normals on the fly
        obj_file_com: Vector3::<f64>::zeros(), // Added later
    };
    // Calculate center of mass, and shift eveything if needed
    let com = calculate_center_of_mass(&rb);
    println!("Center of mass: {:?}", com);
    rb.obj_file_com = com;
    for v in &mut rb.vertices {
        *v -= com;
    }
    // Create rigid particles
    create_rigid_particles(&mut rb, Some(RIGID_BODY_PARTICLES_PER_FACE));
    println!("Rigid particles: {}", rb.rigid_particle_positions.len());

    // Calculate moment of inertia
    let inertia_tensor = calculate_inertia_tensor(&rb);
    rb.inertia_tensor = inertia_tensor;

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
    // let mut vertex_normals = Vec::new();
    // for ind in 0..mesh.normal_indices.len() / 3 {
    //     let v = Vector3::<f64>::new(
    //         mesh.normals[3 * ind] as f64,
    //         mesh.normals[3 * ind + 1] as f64,
    //         mesh.normals[3 * ind + 2] as f64,
    //     );
    //     vertex_normals.push(v);
    // }
    // vertex_normals
    Vec::new()
}

fn create_rigid_particles(rb: &mut RigidBody, num_particles_per_face: Option<usize>) {
    let mut num_p = 0;
    match num_particles_per_face {
        Some(n) => num_p = n,
        None => num_p = RIGID_BODY_PARTICLES_PER_FACE,
    }
    for face_ind in 0..rb.faces.len() {
        let face = rb.faces[face_ind];
        let v1 = rb.vertices[face.0];
        let v2 = rb.vertices[face.1];
        let v3 = rb.vertices[face.2];
        
        for _ in 0..num_p {
            let p = generate_random_point_on_triangle(v1, v2, v3);
            rb.rigid_particle_positions.push(p);
            rb.rigid_particle_triangles.push(face_ind);
            // TODO gonna just use surface normals for now...
            // If triangles are wound incorrectly, this might be wrong
            let face_normal = (v2 - v1).cross(&(v3 - v1)).normalize();
            rb.rigid_particle_normals.push(face_normal);
        }
    }
}