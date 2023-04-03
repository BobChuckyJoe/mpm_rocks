import itertools
import json
import os

import numpy as np
import pandas as pd
import polyscope as ps

TIMESTEP = 0
RENDER_DIST_FIELD = False
ps.init()

# Load the stuff
with open(os.path.expanduser("~") + "/Downloads/sim.json") as f:
    sim = json.load(f)
grid_length = sim["grid_length"]
grid_spacing = sim["grid_spacing"]
num_particles = sim["num_particles"]
particle_positions = sim["particle_positions"]
points = []

# Particles positions
particle_pos = []
for i in range(num_particles):
    particle_pos.append(np.array(particle_positions[TIMESTEP][i]))
particle_point_cloud = ps.register_point_cloud("particle_positions", np.array(particle_pos))

# Generate node points
for i in range(grid_length):
    for j in range(grid_length):
        for k in range(grid_length):
            points.append(np.array([i * grid_spacing, j * grid_spacing, k * grid_spacing]))
points = np.array(points)

ps_cloud = ps.register_point_cloud("grid_nodes", points)

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
