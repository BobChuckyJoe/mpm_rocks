import numpy as np
from scipy.spatial.transform import Rotation as R

def convert_position_to_world(point, orientation, world_position):
    r = R.from_quat(orientation)
    return r.apply(point) + world_position

def convert_vector_to_world(vector, orientation):
    r = R.from_quat(orientation)
    return r.apply(vector)

def get_max_norm(vector_list):
    max_norm = 0
    for vector in vector_list:
        norm = np.linalg.norm(vector)
        if norm > max_norm:
            max_norm = norm
    return max_norm