import json
import os

import numpy as np
import pandas as pd
import polyscope as ps


ps.init()

# Load the stuff
with open(os.path.expanduser("~") + "/Downloads/sim.json") as f:
    sim = json.load(f)
grid_length = sim["grid_length"]
grid_spacing = sim["grid_spacing"]
points = []
vals = []

assert len(sim["unsigned_distance_field"][0]) == grid_length
assert len(sim["unsigned_distance_field"][0][0]) == grid_length
assert len(sim["unsigned_distance_field"][0][0][0]) == grid_length
print(3 * grid_spacing)

# Generate node points
for i in range(grid_length):
    for j in range(grid_length):
        for k in range(grid_length):
           points.append(np.array([i * grid_spacing, j * grid_spacing, k * grid_spacing]))
           vals.append(min(3 * grid_spacing, sim["unsigned_distance_field"][0][i][j][k])) 
points = np.array(points)
vals = np.array(vals)
print(pd.DataFrame(vals).describe())

ps_cloud = ps.register_point_cloud("grid_nodes", points)

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
# ps_cloud.add_scalar_quantity("dist")
ps.show()
