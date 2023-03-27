import bpy
import json
import os
from pathlib import Path
home = str(Path.home())
import time
from mathutils import Quaternion

#with open(os.path.join(home, "Downloads/sim.json")) as f:
#    stuff = json.load(f)

with open("/home/bobchuckyjoe/Projects/mpm_rocks/rigid_body/output.txt") as f:
    lines = f.readlines()

obj = bpy.data.objects["Cylinder"]
obj.rotation_mode = "QUATERNION"
for frame_num, l in enumerate(lines):
    w,x,y,z = map(float, l.split())
    obj.rotation_quaternion = Quaternion((w,x,y,z))
    obj.keyframe_insert(data_path="rotation_quaternion", frame=frame_num)

#box_size = stuff["box_size"]
#grid_length = stuff["grid_length"]
#grid_spacing = stuff["grid_spacing"]
#delta_t = stuff["delta_t"]
#num_particles = stuff["num_particles"]
#num_iterations = stuff["num_iterations"]
#particle_positions = stuff["particle_positions"]

# Select all the spheres and delete them from earlier runs
#for o in bpy.context.scene.objects:
#    if o.type == 'MESH':
#        o.delete()

# Call the operator only once
# bpy.ops.object.delete()

# Each represents a particle
#for _ in range(num_particles):
#    mesh = bpy.ops.mesh.primitive_uv_sphere_add(radius=1, enter_editmode=False, location=(0,0,0), scale=(0.05,0.05,0.05))
        

#print(f"Num iterations: {num_iterations}")
#print(f"Num particles: {num_particles}")
## Do the animation by setting the particle positions
#for t in range(num_iterations):
#    # Skip frames such that it's only 60 fps
#    # TODO

#    curr_time_pos = particle_positions[t]
#    for i in range(num_particles):
#        if i == 0:
#            s = bpy.data.objects["Sphere"]
#        else:
#            s = bpy.data.objects[f"Sphere.{i:03}"]
#        x, y, z = curr_time_pos[i]
##        print(f"{x} {y}")
#        s.location.x = x
#        s.location.y = y
#        s.location.z = z
#        
#        s.keyframe_insert(data_path="location", frame=t)
#    
print("Done!")