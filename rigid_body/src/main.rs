use std::{fs::File, io::Write};

use nalgebra::{Matrix3, Quaternion, UnitQuaternion, Vector3};
use tobj::{self, Mesh};
use tqdm::tqdm;

const DELTA_T: f64 = 0.001;
fn main() {
    let obj_file = "t3.obj";
    let (mut models, materials) =
        tobj::load_obj(obj_file, &tobj::LoadOptions::default()).expect("Failed to load OBJ file");

    println!("Loaded {} models", models.len());
    let m = &mut models[0];
    println!("Model name: {}", m.name);
    let mut mesh = &mut m.mesh;
    println!("Num faces: {}", mesh.face_arities.len());
    println!("num positions: {}", mesh.positions.len() / 3);
    assert!(mesh.positions.len() % 3 == 0);
    println!("face_arities: {:?}", mesh.face_arities);
    println!("Mesh has {} indices", mesh.indices.len());
    // I'm gonna assume these are the indices...
    // println!("Mesh indices: {:?}", mesh.indices);
    println!("Num normals: {}", mesh.normal_indices.len() / 3);
    // Calculate the vertices average
    let mut tot = Vector3::<f64>::zeros();
    for ind in 0..mesh.positions.len() / 3 {
        let v = Vector3::<f64>::new(
            mesh.positions[3 * ind] as f64,
            mesh.positions[3 * ind + 1] as f64,
            mesh.positions[3 * ind + 2] as f64,
        );
        tot += v;
    }
    tot = tot / ((mesh.positions.len() / 3) as f64);
    println!("Average vertex: {:?}", tot);

    // Initial
    let mut omega = Vector3::<f64>::new(0.01, 1.0, 0.01);
    let mut orientation = UnitQuaternion::<f64>::identity();
    let mut position = Vector3::<f64>::zeros(); // Gonna mostly just test rotation around the center of mass

    // Calculate center of mass
    let mut center_of_mass = calculate_center_of_mass(&mesh);
    // Subtract every vertex by the center of mass to center it
    for ind in 0..mesh.positions.len() / 3 {
        mesh.positions[3 * ind] -= center_of_mass.x as f32;
        mesh.positions[3 * ind + 1] -= center_of_mass.y as f32;
        mesh.positions[3 * ind + 2] -= center_of_mass.z as f32;
    }
    println!("Center of mass for mesh: {:?}", center_of_mass);
    // Sanity check
    println!("New com: {:?}", calculate_center_of_mass(&mesh));

    // Calculate inertia tensor
    let inertia_tensor = calculate_inertia_tensor(&mesh);
    println!("Inertia tensor: {:?}", inertia_tensor);
    let svd = inertia_tensor.svd(true, true);
    println!("u: {:?}", svd.u);
    println!("Inertia tensor eigenvalues: {:?}", inertia_tensor.svd(false, false));

    let mut f = File::create("output.txt").unwrap(); // Quaternion at each timestep
    println!("initial L: {:?}", orientation.to_rotation_matrix() * inertia_tensor * orientation.to_rotation_matrix().transpose() * omega);
    for iteration_num in tqdm(0..100000) {
        // let r = orientation.to_rotation_matrix();

        // let r_dot = calculate_star_mat(omega * DELTA_T) * orientation.to_rotation_matrix();
        // TODO: Should it be add? Or should it be multiply?
        // let new_r = r_dot * r.matrix();
        // Update
        // orientation = UnitQuaternion::from_matrix(&new_r);
        let old_orientation = orientation;
        let new_orientation = calculate_new_quat(old_orientation, omega);
        let it = new_orientation.to_rotation_matrix() * inertia_tensor * new_orientation.to_rotation_matrix().inverse();
        
        omega = it.try_inverse().unwrap() * (old_orientation.to_rotation_matrix() * inertia_tensor * old_orientation.to_rotation_matrix().transpose()) * omega;
        orientation = new_orientation;
        if iteration_num % 10 == 0 {
            f.write_all(
                format!(
                    "{} {} {} {}",
                    orientation.w, orientation.i, orientation.j, orientation.k
                )
                .as_bytes(),
            )
            .unwrap();
            f.write_all("\n".as_bytes()).unwrap();
        }
    }
    println!("Final L: {:?}", orientation.to_rotation_matrix() * inertia_tensor * orientation.to_rotation_matrix().transpose() * omega);
}

pub fn calculate_inertia_tensor(mesh: &Mesh) -> Matrix3<f64> {
    let mut inertia_tensor = Matrix3::<f64>::zeros();
    let mut m = 0.0;
    let (mut Cx, mut Cy, mut Cz) = (0.0, 0.0, 0.0); // Centroid
    let (mut xx, mut yy, mut zz, mut yx, mut zx, mut zy) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for triangle_inds in 0..mesh.indices.len() / 3 {
        let vert_a_ind = mesh.indices[3 * triangle_inds] as usize;
        let vert_b_ind = mesh.indices[3 * triangle_inds + 1] as usize;
        let vert_c_ind = mesh.indices[3 * triangle_inds + 2] as usize;
        // let mut vert_inds = (vert_a_ind, vert_b_ind, vert_c_ind);

        let vert_a = Vector3::<f64>::new(
            mesh.positions[3 * vert_a_ind] as f64,
            mesh.positions[3 * vert_a_ind + 1] as f64,
            mesh.positions[3 * vert_a_ind + 2] as f64,
        );
        let vert_b = Vector3::<f64>::new(
            mesh.positions[3 * vert_b_ind] as f64,
            mesh.positions[3 * vert_b_ind + 1] as f64,
            mesh.positions[3 * vert_b_ind + 2] as f64,
        );
        let vert_c = Vector3::<f64>::new(
            mesh.positions[3 * vert_c_ind] as f64,
            mesh.positions[3 * vert_c_ind + 1] as f64,
            mesh.positions[3 * vert_c_ind + 2] as f64,
        );
        // let vertices = [vert_a, vert_b, vert_c];
        // Check the the triangle is oriented correctly
        // if Matrix3::<f64>::from_columns(&[
        //     vert_a,
        //     vert_b,
        //     vert_c,
        // ]).determinant() < 0.0  {
        //     // Swap two vertices
        //     let temp = vert_inds.1;
        //     vert_inds.1 = vert_inds.2;
        //     vert_inds.2 = temp;
        // }

        // let vert_a = vertices[vert_inds.0];
        // let vert_b = vertices[vert_inds.1];
        // let vert_c = vertices[vert_inds.2];
        // Signed tetrahedron volume
        let tet_vol = vert_a.cross(&vert_b).dot(&vert_c); //TODO: I think the code doesn't divide by 6.0 because it's not needed for the inertia tensor...?

        // Contribution to mass
        m += tet_vol;

        // Contribution to centroid
        let sum = vert_a + vert_b + vert_c;
        Cx += tet_vol * sum.x;
        Cy += tet_vol * sum.y;
        Cz += tet_vol * sum.z;

        // Contribution to inertia tensor
        xx += tet_vol
            * (vert_a.x * vert_a.x + vert_b.x * vert_b.x + vert_c.x * vert_c.x + sum.x * sum.x);
        yy += tet_vol
            * (vert_a.y * vert_a.y + vert_b.y * vert_b.y + vert_c.y * vert_c.y + sum.y * sum.y);
        zz += tet_vol
            * (vert_a.z * vert_a.z + vert_b.z * vert_b.z + vert_c.z * vert_c.z + sum.z * sum.z);
        yx += tet_vol
            * (vert_a.y * vert_a.x + vert_b.y * vert_b.x + vert_c.y * vert_c.x + sum.y * sum.x);
        zx += tet_vol
            * (vert_a.z * vert_a.x + vert_b.z * vert_b.x + vert_c.z * vert_c.x + sum.z * sum.x);
        zy += tet_vol
            * (vert_a.z * vert_a.y + vert_b.z * vert_b.y + vert_c.z * vert_c.y + sum.z * sum.y);
    }

    // GetResults function
    let r = 1.0 / (4.0 * m);
    let (Cx, Cy, Cz) = (Cx * r, Cy * r, Cz * r);
    let m = m / 6.0;

    let r = 1.0 / 120.0;
    let Iyx = yx * r - m * Cy * Cx;
    let Izx = zx * r - m * Cz * Cx;
    let Izy = zy * r - m * Cz * Cy;

    let xx = xx * r - m * Cx * Cx;
    let yy = yy * r - m * Cy * Cy;
    let zz = zz * r - m * Cz * Cz;

    let Ixx = yy + zz;
    let Iyy = zz + xx;
    let Izz = xx + yy;

    inertia_tensor[(0, 0)] = Ixx;
    inertia_tensor[(1, 1)] = Iyy;
    inertia_tensor[(2, 2)] = Izz;
    inertia_tensor[(0, 1)] = Iyx;
    inertia_tensor[(1, 0)] = Iyx;
    inertia_tensor[(0, 2)] = Izx;
    inertia_tensor[(2, 0)] = Izx;
    inertia_tensor[(1, 2)] = Izy;
    inertia_tensor[(2, 1)] = Izy;

    inertia_tensor
}

fn calculate_center_of_mass(m: &Mesh) -> Vector3<f64> {
    let mut tot_vol = 0.0;
    let mut center_sum = Vector3::<f64>::new(0.0, 0.0, 0.0);
    for ind_start in 0..m.indices.len() / 3 {
        let ind_a = m.indices[3 * ind_start] as usize;
        let ind_b = m.indices[3 * ind_start + 1] as usize;
        let ind_c = m.indices[3 * ind_start + 2] as usize;
        let vert_a = Vector3::<f64>::new(
            m.positions[3 * ind_a] as f64,
            m.positions[3 * ind_a + 1] as f64,
            m.positions[3 * ind_a + 2] as f64,
        );
        let vert_b = Vector3::<f64>::new(
            m.positions[3 * ind_b] as f64,
            m.positions[3 * ind_b + 1] as f64,
            m.positions[3 * ind_b + 2] as f64,
        );
        let vert_c = Vector3::<f64>::new(
            m.positions[3 * ind_c] as f64,
            m.positions[3 * ind_c + 1] as f64,
            m.positions[3 * ind_c + 2] as f64,
        );

        let center = (vert_a + vert_b + vert_c) / 4.0; // The origin is implicitly the 4th poitn
        let volume = vert_a.dot(&vert_b.cross(&vert_c)) / 6.0;
        tot_vol += volume;
        center_sum += center * volume;
    }
    center_sum / tot_vol
}

fn calculate_star_mat(v: Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0)
}

fn calculate_new_quat(old_q: UnitQuaternion<f64>, omega: Vector3<f64>) -> UnitQuaternion<f64> {
    let omega = omega * DELTA_T;
    let rest = omega / omega.norm() * (omega.norm() / 2.0).sin();
    let q = UnitQuaternion::from_quaternion(Quaternion::new(
        (0.5 * omega.norm()).cos(),
        rest.x,
        rest.y,
        rest.z,
    ));
    q * old_q
}
