use std::fs::File;
use std::io::{BufWriter, Write};

use nalgebra::Vector3;

use crate::types::{RigidBody, Particle};
/// Only output rigid body and particle data
pub struct FileOutput {
    output_directory: String,
    frame_counter: usize,
    original_obj_file_string: Vec<String>,
}

impl FileOutput {
    pub fn new(output_directory: String, rigid_body_obj_txt: String) -> FileOutput {
        FileOutput {
            output_directory,
            frame_counter: 0,
            original_obj_file_string: rigid_body_obj_txt.split("\n").map(|s| s.to_string()).collect::<Vec<_>>(),
        }
    }
    pub fn write_frame(&mut self, rb: &RigidBody, particles: &Vec<Particle>) {
        // Write particle data
        let mut particle_file = File::create(format!("{}/{}.csv", self.output_directory, self.frame_counter)).unwrap();
        for p in particles.iter() {
            writeln!(particle_file, "{},{},{}", p.position.x, p.position.y, p.position.z).unwrap();
        }
        // Write rigid body position and orientations to obj files
        let mut rigid_body_file = BufWriter::new(File::create(format!("{}/{}.obj", self.output_directory, self.frame_counter)).unwrap());
        
        for l in &self.original_obj_file_string {
            let line_split = l.split(" ").collect::<Vec<_>>();
            if l.starts_with("vn") && line_split.len() == 4 {
                let vn = Vector3::<f64>::new(line_split[1].parse::<f64>().unwrap(),
                                              line_split[2].parse::<f64>().unwrap(),
                                              line_split[3].parse::<f64>().unwrap());
                let vn = rb.orientation * vn;
                writeln!(rigid_body_file, "vn {} {} {}", vn.x, vn.y, vn.z).unwrap();
            }
            else if l.starts_with("v") && line_split.len() == 4 {
                let v = Vector3::<f64>::new(line_split[1].parse::<f64>().unwrap(),
                                             line_split[2].parse::<f64>().unwrap(),
                                             line_split[3].parse::<f64>().unwrap());
                let v = rb.orientation * v + rb.position;
                writeln!(rigid_body_file, "v {} {} {}", v.x, v.y, v.z).unwrap();
            }
            else {
                writeln!(rigid_body_file, "{}", l).unwrap();
            }
        }
        self.frame_counter += 1;
    }
}