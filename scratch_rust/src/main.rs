use rand::Rng;
use nalgebra::{Matrix3};
use tqdm::tqdm;

fn generate_random_mat() -> Matrix3<f64> {
    let mut rng = rand::thread_rng();
    let mut mat = Matrix3::zeros();
    for i in 0..3 {
        for j in 0..3 {
            mat[(i, j)] = rng.gen();
        }
    }
    mat
}
fn get_singular_matrix(a: f64, b:f64, c:f64) -> Matrix3<f64> {
    let mut mat = Matrix3::zeros();
    mat[(0, 0)] = a;
    mat[(1,1)] = b;
    mat[(2,2)] = c;
    mat
}

fn main() {
    for i in tqdm(0..10000000) {
        let mat = generate_random_mat();
        // println!("mat: {}", mat);
        let svd = mat.svd(true, true);
        assert!((svd.recompose().unwrap() - mat).norm() < 1e-6);
        let singular_values = get_singular_matrix(svd.singular_values[0], svd.singular_values[1], svd.singular_values[2]);
        assert!((svd.u.unwrap() * singular_values * svd.v_t.unwrap() - mat).norm() < 1e-6);
    }
}
