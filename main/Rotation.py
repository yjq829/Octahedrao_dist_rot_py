# ******************************************************************************************************
# python 3.6
# version 1.0
#
# object of octahedra part: *.vasp file
# (1) find 8 different octahedra
# (2) define each octahedra
# (3) all data restored is atom index ID in atoms
# ['C', 'N', 'H', 'Ge', 'Cl', 'Br', 'I', 'Sn', 'Pb', 'Cs', 'Ba', 'Sr', 'Ca', 'Rb', 'K']
#
# ******************************************************************************************************


# dependency environment
import math

import numpy as np


# calculate rotation information from oct_ref and oct_test
# rotation degrees: compare angel between vector(x_ver)
# input - oct_ref, oct_test
# output - deg_rot
def rotation_degree(oct_ref, oct_test):
    vec_ref = oct_ref.x_ver_positions_c[0] - oct_ref.x_ver_positions_c[1]
    vec_ref_norm = vec_ref / np.linalg.norm(vec_ref)
    vec_test = oct_test.x_ver_positions_c[0] - oct_test.x_ver_positions_c[1]
    vec_test_norm = vec_test / np.linalg.norm(vec_test)
    rot_ang = np.arccos(np.dot(vec_ref_norm, vec_test_norm))
    rot_ang = abs(rot_ang * 180 / math.pi)
    if rot_ang > 90:
        rot_ang = 180 - rot_ang
    return rot_ang
