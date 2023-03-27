
use std::collections::HashSet;

use nalgebra::{Matrix3, UnitQuaternion, Vector3};

use crate::config::RIGID_BODY_PARTICLES_PER_FACE;
use crate::equations::{calculate_inertia_tensor, generate_random_point_on_triangle};
use crate::types::RigidBody;

fn create_vec3(x: f64, y: f64, z: f64) -> Vector3<f64> {
    Vector3::<f64>::new(x, y, z)
}



// TODO (rust thing) learn how to make this a constant
fn get_vertices() -> Vec<Vector3<f64>> {
    vec![
        create_vec3(0.019649, 0.840808, 0.448559),
        create_vec3(0.828664, 0.531788, -0.051401),
        create_vec3(0.019649, 0.840809, -0.551401),
        create_vec3(-0.480351, 0.031790, 0.757579),
        create_vec3(0.019649, -0.777226, 0.448559),
        create_vec3(0.828664, -0.468209, -0.051401),
        create_vec3(-0.789368, 0.531788, -0.051401),
        create_vec3(-0.480351, 0.031791, -0.860501),
        create_vec3(0.019649, -0.777226, -0.551401),
        create_vec3(-0.789368, -0.468209, -0.051401),
        create_vec3(0.519654, 0.031790, 0.757579),
        create_vec3(0.519654, 0.031791, -0.860501),
    ]
}

fn get_center_of_mass() -> Vector3<f64> {
    let mut sum = Vector3::<f64>::zeros();
    let vertices = get_vertices();
    for vertex in vertices.iter() {
        sum += vertex;
    }
    let center_of_mass = sum / vertices.len() as f64;
    center_of_mass
}

pub fn create_icosahedron() -> RigidBody {
    let center_of_mass = get_center_of_mass();
    println!("center of mass: {:?}", center_of_mass);
    // Sanity check that the center of mass is correct (i.e. equidistant to all of the vertices)
    for vertex in get_vertices().iter() {
        println!("distance: {}", (vertex - center_of_mass).norm());
    }
    // ^ Looks good to me

    // Get all of the vertices to be in the "material frame" (with the origin at the center of mass)
    let mut vertices: Vec<Vector3<f64>> = get_vertices();
    for vertex in vertices.iter_mut() {
        *vertex -= center_of_mass;
    }
    // Generate all of the faces
    let mut vert_faces: HashSet<(usize, usize, usize)> = HashSet::new();
    for v_ind in 0..vertices.len() {
        let mut other_inds: Vec<(usize, f64)> = vec![];
        for other_v_ind in 0..vertices.len() {
            if v_ind == other_v_ind {
                continue;
            }
            let distance = (vertices[v_ind] - vertices[other_v_ind]).norm();
            other_inds.push((other_v_ind, distance));
        }
        // Sort by distance
        other_inds.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        assert!(other_inds[0].1 < other_inds[other_inds.len() - 1].1);
        // Since each vertex is connected to 5 other vertices, we can just take the first 5
        for i in 0..5 {
            for j in 0..5 {
                if i == j {
                    continue;
                }
                assert!((other_inds[i].1 - other_inds[j].1).abs() < 1e-3);
                // Check if they belong in the same triangle
                if ((vertices[v_ind] - vertices[other_inds[i].0]).norm()
                    - (vertices[other_inds[i].0] - vertices[other_inds[j].0]).norm())
                .abs()
                    > 1e-3
                {
                    continue;
                }
                if (vertices[v_ind] - vertices[other_inds[j].0]).norm()
                    - (vertices[other_inds[i].0] - vertices[other_inds[j].0]).norm()
                    > 1e-3
                {
                    continue;
                }
                let mut face_verts = vec![v_ind, other_inds[i].0, other_inds[j].0];
                face_verts.sort(); // Sort by ascending index
                vert_faces.insert((face_verts[0], face_verts[1], face_verts[2]));
            }
        }
    }

    for (a, b, c) in vert_faces.iter() {
        let p1 = vertices[*a];
        let p2 = vertices[*b];
        let p3 = vertices[*c];
        assert!(
            ((p1 - p2).norm() - (p2 - p3).norm()).abs() < 1e-3,
            "p1 - p2: {:?}, p2 - p3: {:?}",
            (p1 - p2).norm(),
            (p2 - p3).norm()
        );
        assert!(
            ((p1 - p2).norm() - (p1 - p3).norm()).abs() < 1e-3,
            "p1 - p2: {:?}, p1 - p3: {:?}",
            (p1 - p2).norm(),
            (p1 - p3).norm()
        );
    }
    println!("vert faces: {:?}", vert_faces);
    assert!(
        vert_faces.len() == 20,
        "vert_faces.len() = {}",
        vert_faces.len()
    );
    // Convert it into a list
    let vert_face_inds: Vec<(usize, usize, usize)> = vert_faces.into_iter().collect();
    let face_normals: Vec<Vector3<f64>> = vert_face_inds
        .iter()
        .map(|(a, b, c)| {
            let p1 = vertices[*a];
            let p2 = vertices[*b];
            let p3 = vertices[*c];
            let normal = (p2 - p1).cross(&(p3 - p1)).normalize();
            normal
        })
        .collect();
    // Sample particles on each face
    // I'm using uniform random, but the paper uses a lattice
    let mut rigid_particle_positions: Vec<Vector3<f64>> = vec![];
    let mut rigid_particle_triangles: Vec<usize> = vec![]; // Indices into the vertices array
    let mut rigid_particle_normals: Vec<Vector3<f64>> = vec![];
    for face_ind in 0..vert_face_inds.len() {
        let (face_ind_1, face_ind_2, face_ind_3) = vert_face_inds[face_ind];
        let p1 = vertices[face_ind_1];
        let p2 = vertices[face_ind_2];
        let p3 = vertices[face_ind_3];
        for _ in 0..RIGID_BODY_PARTICLES_PER_FACE {
            let random_point = generate_random_point_on_triangle(p1, p2, p3);
            rigid_particle_positions.push(random_point);
            rigid_particle_triangles.push(face_ind);
            // Calculate surface normal
            let triangle_centroid: Vector3<f64> = (p1 + p2 + p3) / 3.0;
            let derp = triangle_centroid - p2;
            let mut poss_normal: Vector3<f64> = (triangle_centroid - p1).cross(&derp);
            if poss_normal.dot(&triangle_centroid) < 0.0 {
                poss_normal *= -1.0;
            }
            rigid_particle_normals.push(poss_normal / poss_normal.norm());
        }
    }
    // TODO I'm just gonna use moment of inertial for a sphere
    let radius = rigid_particle_positions[0].norm();
    let mass = 1.0;

    let mut rb = RigidBody {
        position: Vector3::<f64>::zeros(),
        velocity: Vector3::<f64>::zeros(),
        orientation: UnitQuaternion::<f64>::identity(),
        omega: Vector3::<f64>::zeros(),
        mass: mass,
        vertices: vertices,
        vertex_normals: vec![], // TODO
        faces: vert_face_inds,
        rigid_particle_positions: rigid_particle_positions,
        rigid_particle_triangles: rigid_particle_triangles,
        rigid_particle_normals: rigid_particle_normals,
        moment_of_inertia: Matrix3::zeros(),
        obj_file_com: Vector3::<f64>::zeros(),
    };
    let inertial_tensor = calculate_inertia_tensor(&rb);
    rb.moment_of_inertia = inertial_tensor;
    rb
}

/*
v 0.019649 0.840808 0.448559
v 0.828664 0.531788 -0.051401
v 0.019649 0.840809 -0.551401
v -0.480351 0.031790 0.757579
v 0.019649 -0.777226 0.448559
v 0.828664 -0.468209 -0.051401
v -0.789368 0.531788 -0.051401
v -0.480351 0.031791 -0.860501
v 0.019649 -0.777226 -0.551401
v -0.789368 -0.468209 -0.051401
v 0.519654 0.031790 0.757579
v 0.519654 0.031791 -0.860501
 */
