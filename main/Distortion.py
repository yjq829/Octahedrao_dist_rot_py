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


# calculate distortion information from oct_ref and oct_test
# angle distortion: compare angels of Xver-B-Xver
# vertical distortion: compare average distance of Xver-B
# horizonal distortion: compare average distace of Xhor-B
# input - oct_ref, oct_test
# output - ang_dist
#           ver_dist
#           hor_dist
def ang_distortion(oct_ref, oct_test):
    # calculate ref angle
    vec_ref_1 = oct_ref.x_ver_positions_c[0] - oct_ref.atoms.positions[oct_ref.b_center]
    vec_ref_1n = vec_ref_1 / np.linalg.norm(vec_ref_1)
    vec_ref_2 = oct_ref.x_ver_positions_c[1] - oct_ref.atoms.positions[oct_ref.b_center]
    vec_ref_2n = vec_ref_2 / np.linalg.norm(vec_ref_2)
    ang_ref = np.arccos(np.dot(vec_ref_1n, vec_ref_2n))
    # print(ang_ref*180/math.pi)
    # calculate test angle
    vec_test_1 = oct_test.x_ver_positions_c[0] - oct_test.atoms.positions[oct_test.b_center]
    vec_test_1n = vec_test_1 / np.linalg.norm(vec_test_1)
    vec_test_2 = oct_test.x_ver_positions_c[1] - oct_test.atoms.positions[oct_test.b_center]
    vec_test_2n = vec_test_2 / np.linalg.norm(vec_test_2)
    ang_test = np.arccos(np.dot(vec_test_1n, vec_test_2n))
    # print(ang_test*180/math.pi)
    dis_ang = abs((ang_test - ang_ref) * 180 / math.pi)
    if dis_ang > 90:
        dis_ang = 180 - dis_ang
    return dis_ang


def ver_dist(oct_ref, oct_test):
    # calculate ver distance average for ref
    vec_ref_1 = np.linalg.norm(oct_ref.x_ver_positions_c[0] - oct_ref.atoms.positions[oct_ref.b_center])
    vec_ref_2 = np.linalg.norm(oct_ref.x_ver_positions_c[1] - oct_ref.atoms.positions[oct_ref.b_center])
    ver_dist_ref = (vec_ref_1 + vec_ref_2) / 2
    # calculate ver distances average for test
    vec_test_1 = np.linalg.norm(oct_test.x_ver_positions_c[0] - oct_test.atoms.positions[oct_test.b_center])
    vec_test_2 = np.linalg.norm(oct_test.x_ver_positions_c[1] - oct_test.atoms.positions[oct_test.b_center])
    ver_dist_test = (vec_test_1 + vec_test_2) / 2
    return ver_dist_test - ver_dist_ref


def hor_dist(oct_ref, oct_test):
    hor_dis_ref = []
    hor_dis_test = []
    for i in range(len(oct_ref.x_hor_positions_c)):
        hor_dis_ref.append(np.linalg.norm(oct_ref.x_hor_positions_c[i] - oct_ref.atoms.positions[oct_ref.b_center]))
    for i in range(len(oct_test.x_hor_positions_c)):
        hor_dis_test.append(np.linalg.norm(oct_test.x_hor_positions_c[i] - oct_test.atoms.positions[oct_test.b_center]))

    hor_avg_ref = sum(hor_dis_ref) / len(hor_dis_ref)
    # print(hor_avg_ref)
    hor_avg_test = sum(hor_dis_test) / len(hor_dis_test)
    # print(hor_avg_test)
    return hor_avg_test - hor_avg_ref
