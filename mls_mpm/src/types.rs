use nalgebra::{Matrix3, UnitQuaternion, Vector3};
#[derive(Copy, Clone)]
pub struct Particle {
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    // Since we're using the cubic interpolation stencils, the D matrix is going to be a constant
    // Therefore, we will only need to store the B matrix.
    // C = B * D^{-1}
    pub apic_b: Matrix3<f64>,
    pub deformation_gradient: Matrix3<f64>,
    pub mass: f64,
    pub density: f64,
    pub tag: i32,                      // This is the T_{pr} from equation 21
    pub particle_distance: f64, // d_p from 5.3.2 of MLS MPM paper. I'm also assuming that sign(d_p) = T_{ir}
    pub particle_normal: Vector3<f64>, // The particle normal from 5.3.2 of MLS MPM paper
}

#[derive(Copy, Clone)]
pub struct Gridcell {
    pub velocity: Vector3<f64>,
    pub mass: f64,
    pub force: Vector3<f64>,
    pub unsigned_distance: f64, // Smallest unsigned distance to a rigid particle
    pub rigid_particle_index: i32, // the index of the rigid particle that this gridcell is closest to
    pub distance_sign: i32,        // T_{ir}
}

impl Gridcell {
    pub fn new() -> Gridcell {
        Gridcell {
            velocity: Vector3::zeros(),
            mass: 0.0,
            force: Vector3::zeros(),
            unsigned_distance: 1000000.0,
            rigid_particle_index: -1,
            distance_sign: 0,
        }
    }

    pub fn reset_values(&mut self) {
        self.velocity.scale_mut(0.0);
        self.mass = 0.0;
        self.force.scale_mut(0.0);
        self.unsigned_distance = 1000000.0; // Because, we're doing a min operation, this should be large at first
        self.distance_sign = 0;
    }
}

pub struct RigidBody {
    // Position and velocity are in world coords
    pub position: Vector3<f64>, // Position of center of mass
    pub velocity: Vector3<f64>,
    // These are in material coords
    pub orientation: UnitQuaternion<f64>, // Rotation matrix. Should I use quaternions?
    pub angular_momentum: Vector3<f64>, // Angular MOMENTUM
    pub mass: f64,
    pub inertia_tensor: Matrix3<f64>, // Also in material coords
    
    // Mesh data
    pub vertices: Vec<Vector3<f64>>,
    pub faces: Vec<(usize, usize, usize)>, // Triangles made up of indices of the vertices
    pub vertex_normals: Vec<Vector3<f64>>, // Per vertex normals
    pub obj_file_com: Vector3<f64>,

    // The particles will be in the material frame
    pub rigid_particle_positions: Vec<Vector3<f64>>, // TODO should i use particles or should i just use positions...?
    pub rigid_particle_triangles: Vec<usize>, // Corresponding triangles for each particle. Indexes into the vertices field
    pub rigid_particle_normals: Vec<Vector3<f64>>, // Corresponding normals for each triangle
}

pub enum BoundaryCondition {
    STICKY,
    SLIP,
    SEPARATE,
}
