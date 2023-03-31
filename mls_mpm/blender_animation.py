import bpy
import itertools
import json
import os
from pathlib import Path
home = str(Path.home())
import time
from mathutils import Quaternion

with open(os.path.join(home, "Downloads/sim.json")) as f:
    stuff = json.load(f)

box_size = stuff["box_size"]
grid_length = stuff["grid_length"]
grid_spacing = stuff["grid_spacing"]
delta_t = stuff["delta_t"]
num_particles = stuff["num_particles"]
num_iterations = stuff["num_iterations"]
particle_positions = stuff["particle_positions"]
# Rigid body stuff
rigid_body_positions = stuff["rigid_body_positions"]
rigid_body_orientations = stuff["rigid_body_orientations"]
rigid_body_velocities = stuff["rigid_body_velocities"]
rigid_body_omegas = stuff["rigid_body_omegas"]
obj_file_com = stuff["obj_file_com"]
unsigned_distance_field = stuff["unsigned_distance_field"]

rigid_body = bpy.data.objects["rigid_body"]

def rgb(minimum, maximum, value):
    minimum, maximum = float(minimum), float(maximum)
    ratio = 2 * (value-minimum) / (maximum - minimum)
    b = int(max(0, 255*(1 - ratio)))
    r = int(max(0, 255*(ratio - 1)))
    g = 255 - b - r
    return r, g, b

def init():
    # Each represents a particle
    for _ in range(num_particles):
        mesh = bpy.ops.mesh.primitive_uv_sphere_add(radius=1, enter_editmode=False, location=(0,0,0), scale=(0.05,0.05,0.05))
    # TODO For now, just manually load the rigid body with the name "rigid_body"
    # Set the origin of the rigid_body
    bpy.context.scene.cursor.location = obj_file_com
    rigid_body.select_set(True)
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR', center="MEDIAN")
    # Set rotation mode to quaternion
    rigid_body.rotation_mode = 'QUATERNION'
    # Add balls to represent unsigned distance field values
    for i,j,k in itertools.product(range(len(unsigned_distance_field)), repeat=3):
        mesh = bpy.ops.mesh.primitive_uv_sphere_add(name=f"grid_unsigned_{i}_{j}_{k}", radius=0.01, enter_editmode=False, location=(i * grid_spacing, j * grid_spacing, k * grid_spacing), scale=(1, 1, 1))
        # bpy.data.objects.new(f"grid_unsigned_{i}_{j}_{k}", mesh)
        # mat = bpy.data.materials.new(name="derp")
        # bpy.data.objects[f"grid_unsigned_{i}_{j}_{k}"].active_material = mat
init()

print(f"Num iterations: {num_iterations}")
print(f"Num particles: {num_particles}")
# Do the animation by setting the particle positions
for t in range(num_iterations):
    # Skip frames such that it's only 60 fps
    # TODO

    curr_time_pos = particle_positions[t]
    for i in range(num_particles):
        if i == 0:
            s = bpy.data.objects["Sphere"]
        else:
            s = bpy.data.objects[f"Sphere.{i:03}"]
        x, y, z = curr_time_pos[i]
#        print(f"{x} {y}")
        s.location.x = x
        s.location.y = y
        s.location.z = z
        
        s.keyframe_insert(data_path="location", frame=t)
    curr_time_rigid_body_pos = rigid_body_positions[t]
    curr_time_rigid_body_orientation = rigid_body_orientations[t]
    rigid_body.location = curr_time_rigid_body_pos
    rigid_body.rotation_quaternion = curr_time_rigid_body_orientation
    rigid_body.keyframe_insert(data_path="location", frame=t)
    rigid_body.keyframe_insert(data_path="rotation_quaternion", frame=t)
    
    # # Change grid unsigned distance field colors
    # for i,j,k in itertools.product(range(len(unsigned_distance_field)), repeat=3):
    #     grid_unsigned = bpy.data.objects[f"grid_unsigned_{i}_{j}_{k}"]
    #     grid_unsigned.color = rgb(0, 3 * grid_spacing, unsigned_distance_field[i][j][k])
    #     grid_unsigned.keyframe_insert(data_path="color", frame=t)
print("Done!")