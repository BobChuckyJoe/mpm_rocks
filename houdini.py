from scipy.spatial.transform import Rotation as R
import numpy as np
from tqdm import tqdm

with open("sim.json") as f:
    import json
    data = json.load(f)
print("Done loading sim")
with open("icosahedron.obj") as f:
    rigid_body_lines = f.readlines()
print("Done loading obj")
def process_obj(lines, rot_mat, position):
    out = ""
    for l in lines:
        if l.startswith("vn") and len(l.split()) == 4:
            split = l.split()
            vn = np.array([float(split[1]), float(split[2]), float(split[3])])
            new_vn = rot_mat @ vn 
            # print(f"rot_mat: {rot_mat}")
            # print(f"vn: {vn}")
            # print(f"new_vn: {new_vn}")
            out += f"vn {new_vn[0]} {new_vn[1]} {new_vn[2]}\n"
        elif l.startswith("v") and len(l.split()) == 4:
            split = l.split()
            v = np.array([float(split[1]), float(split[2]), float(split[3])])
            # print(rot_mat)
            # print(v)
            # print(position)
            new_v = rot_mat @ v + position
            out += f"v {new_v[0]} {new_v[1]} {new_v[2]}\n"
        else:
            out += l
    return out
for t in tqdm(range(len(data["particle_positions"]))):
    # Output particle positions
    with open(f"{t}.csv", "w") as f:
        for p in data["particle_positions"][t]:
            f.write(f"{p[0]},{p[1]},{p[2]}\n")
    # Output rigid_body OBJs
    with open(f"{t}.obj", "w") as f:
        position, orientation = (data["rigid_body_positions"][t], data["rigid_body_orientations"][t])
        rot = R.from_quat(orientation).as_matrix()
        position = np.array(position)
        f.write(process_obj(rigid_body_lines, rot, position))
        