import itertools
import json
import os

import numpy as np
import pandas as pd
import polyscope as ps

from math_utils import *

TIMESTEP = 0
RENDER_DIST_FIELD = False
ps.init()
# Consistent with blender
ps.set_up_dir("z_up")

# Load the stuff
with open(os.path.expanduser("~") + "/Downloads/sim.json") as f:
    sim = json.load(f)

grid_affinities = sim["grid_affinities"]
grid_distance_signs = sim["grid_distance_signs"]
grid_length = sim["grid_length"]
grid_spacing = sim["grid_spacing"]
num_particles = sim["num_particles"]
particle_positions = sim["particle_positions"]
rigid_body_positions = sim["rigid_body_positions"]
curr_rigid_body_position = rigid_body_positions[TIMESTEP]
rigid_body_orientations = sim["rigid_body_orientations"]
curr_rigid_body_orientation = rigid_body_orientations[TIMESTEP]
rigid_body_velocities = sim["rigid_body_velocities"]
rigid_body_angular_momentums = sim["rigid_body_angular_momentums"]
rigid_body_vertices = sim["rigid_body_vertices"]
rigid_body_triangles = sim["rigid_body_triangles"]
rigid_body_coms = sim["obj_file_com"]
rigid_particle_positions = sim["rigid_particle_positions"]
rigid_particle_triangles = sim["rigid_particle_triangles"]

# Particles positions
particle_pos = []
for i in range(num_particles):
    particle_pos.append(np.array(particle_positions[TIMESTEP][i]))
particle_point_cloud = ps.register_point_cloud("particle_positions", np.array(particle_pos))

# Generate node points
points = []
for i in range(grid_length):
    for j in range(grid_length):
        for k in range(grid_length):
            points.append(np.array([i * grid_spacing, j * grid_spacing, k * grid_spacing]))
points = np.array(points)
ps_cloud = ps.register_point_cloud("grid_nodes", points)

# Rigid body mesh
world_vert_pos = []
for i in range(len(rigid_body_vertices)):
    world_vert_pos.append(np.array(convert_position_to_world(rigid_body_vertices[i], curr_rigid_body_orientation, curr_rigid_body_position)))
rigid_body_mesh = ps.register_surface_mesh("rigid_body", np.array(world_vert_pos), np.array(rigid_body_triangles))

# Rigid body particles
rigid_particle_position_points = []
for i in range(len(rigid_particle_positions)):
    rigid_particle_position_points.append(np.array(convert_position_to_world(rigid_particle_positions[i], curr_rigid_body_orientation, curr_rigid_body_position)))
rigid_particle_position_points = np.array(rigid_particle_position_points)
rigid_particle_position_point_cloud = ps.register_point_cloud("rigid_particle_positions", rigid_particle_position_points)

# Rigid body particle face normals (have vectors that point out of each point)
rigid_particle_face_normals = []
for i in range(len(rigid_particle_triangles)):
    va, vb, vc = rigid_body_triangles[rigid_particle_triangles[i]]
    va = np.array(rigid_body_vertices[va])
    vb = np.array(rigid_body_vertices[vb])
    vc = np.array(rigid_body_vertices[vc])
    ab = vb - va
    ac = vc - va
    normal = np.cross(ab, ac)
    normal = normal / np.linalg.norm(normal)
    rigid_particle_face_normals.append(normal)
rigid_particle_face_normals = np.array(rigid_particle_face_normals)
rigid_particle_position_point_cloud.add_vector_quantity("rigid_particle_face_normals", rigid_particle_face_normals, radius = 0.001, length = 0.005)

# Grid affinities, which shows which nodes have a valid distance to the rigid body
grid_affinities_locs = []
for i in range(grid_length):
    for j in range(grid_length):
        for k in range(grid_length):
            if grid_affinities[TIMESTEP][i][j][k]:
                grid_affinities_locs.append(np.array([i * grid_spacing, j * grid_spacing, k * grid_spacing]))
ps.register_point_cloud("grid_affinities", np.array(grid_affinities_locs))

# Distance signs
neg_dist = []
for i in range(grid_length):
    for j in range(grid_length):
        for k in range(grid_length):
            if grid_distance_signs[TIMESTEP][i][j][k] == -1:
                neg_dist.append(np.array([i * grid_spacing, j * grid_spacing, k * grid_spacing]))
ps.register_point_cloud("neg_dist", np.array(neg_dist))

pos_dist = []
for i in range(grid_length):
    for j in range(grid_length):
        for k in range(grid_length):
            if grid_distance_signs[TIMESTEP][i][j][k] == 1:
                pos_dist.append(np.array([i * grid_spacing, j * grid_spacing, k * grid_spacing]))
ps.register_point_cloud("pos_dist", np.array(pos_dist))

zero_dist = []
for i in range(grid_length):
    for j in range(grid_length):
        for k in range(grid_length):
            if grid_distance_signs[TIMESTEP][i][j][k] == 0:
                zero_dist.append(np.array([i * grid_spacing, j * grid_spacing, k * grid_spacing]))
ps.register_point_cloud("zero_dist", np.array(zero_dist))

if RENDER_DIST_FIELD:
    vals = []
    for i,j,k in itertools.product(range(grid_length), repeat=3):
        vals.append(min(3 * grid_spacing, sim["unsigned_distance_field"][TIMESTEP][i][j][k])) 
    vals = np.array(vals)
    print(pd.DataFrame(vals).describe())
    ps_cloud.add_scalar_quantity("dist", vals, vminmax=(0, 0.06), cmap="viridis")
    
    # Set radius to 0 for points that have no distance
    rads = []
    for i in vals:
        if i == 3 * grid_spacing:
            rads.append(0)
        else:
            rads.append(0.0001)
        # rads.append(0.0001)
    rads= np.array(rads)

    ps_cloud.add_scalar_quantity("radii", rads)
    ps_cloud.set_point_radius_quantity("radii")

# Grid velocities
grid_vels = []
for i in range(grid_length):
    for j in range(grid_length):
        for k in range(grid_length):
            grid_v = np.array(sim["grid_velocities"][TIMESTEP][i][j][k])
            if grid_v[0] != 0 or grid_v[1] != 0 or grid_v[2] != 0:
                print(f"{grid_v} at {i}, {j}, {k}")
            grid_vels.append(np.array(sim["grid_velocities"][TIMESTEP][i][j][k]) * 1e5)
grid_vels = np.array(grid_vels)
ps_cloud.add_vector_quantity("grid_vels", grid_vels, radius = 0.001, length = 0.005)

# ps_cloud.add_scalar_quantity("dist")
ps.show()
