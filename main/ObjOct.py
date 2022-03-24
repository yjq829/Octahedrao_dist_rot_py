# ******************************************************************************************************
# python 3.6
# version 1.0
#
# object of octahedral part: *.vasp file
# (1) find 8 different octahedral
# (2) define each octahedral
# (3) all data restored is atom index ID in atoms
# ['C', 'N', 'H', 'Ge', 'Cl', 'Br', 'I', 'Sn', 'Pb', 'Cs', 'Ba', 'Sr', 'Ca', 'Rb', 'K']
#
# ******************************************************************************************************

# dependency environment
from ase.io import read
import ase.neighborlist
import numpy as np
import heapq


# oct_obj:
# b_center - center B atom index
# x_hor - a list of 4 X atoms that on the horizonal plane
# x_ver - a list of 2 X atoms that in vertical direction
# x_*_positions_c - a list of positions (ndarry) corrected by periodic boundary
# atoms - the atoms obj from reading the files
class OctObj:
    def __init__(self, b_center, x_hor_list, x_ver_list, x_hor_positions_c, x_ver_positions_c, atoms):
        self.b_center = b_center
        self.x_hor = x_hor_list
        self.x_ver = x_ver_list
        self.atoms = atoms
        self.x_hor_positions_c = x_hor_positions_c
        self.x_ver_positions_c = x_ver_positions_c

    def atom_ver_distance(self):
        a1 = np.array(self.atoms.get_positions()[x_ver[0]])
        a2 = np.array(self.atoms.get_positions()[x_ver[1]])
        return np.linalg.norm(a1 - a2)


# rugular used fundions


# B atoms finding funciton
# input - atoms obj read from POSCAR
# output - a index list of B atoms
def b_finding(atoms):
    ind_ele = enumerate(atoms.get_chemical_symbols())
    b_index = [i for i, x in ind_ele if x in ["Ge", "Pb", "Sn", "Ba", "Sr", "Ca"]]
    return b_index


# X atoms finding funciton
# input - atoms obj read from POSCAR
# output - a index list of x atoms
def x_finding(atoms):
    ind_ele = enumerate(atoms.get_chemical_symbols())
    x_index = [i for i, x in ind_ele if x in ["Cl", "Br", "I"]]
    return x_index


# X atoms in vertical finding function
# input - center B atom id, x atom id, atoms obj
# output - 6 neighbor x atoms with first and second min distance with given B atom

def find_x_oct(b_center_id, atoms):
    # using ase neighbor module
    n_cut = ase.neighborlist.natural_cutoffs(atoms, mult=1.5)
    # to get the secondary neighbor, increase cutoff by 0.5A
    # n_cut=[n_cut[i]+0.4 for i in range(len(n_cut))]
    n1 = ase.neighborlist.NeighborList(n_cut, self_interaction=False, bothways=True)
    n1.update(atoms)
    x_oct_finding, x_offsets = n1.get_neighbors(b_center_id)
    x_oct_finding = list(x_oct_finding)

    # checking if A atoms is included in x_oct_finding
    x_index_using = x_finding(atoms)
    x_oct_real = []
    x_offsets_real = []
    for x in x_oct_finding:
        if x in x_index_using:
            x_oct_real.append(x)
            x_offsets_real.append(x_offsets[x_oct_finding.index(x)])
    # if b_center_id in x_oct:
    #     x_oct.remove(b_center_id)
    return x_oct_real, x_offsets_real


# finding two types of X atoms, horizonal x atoms are closer to average z coordination
# input - center B id, x id in octahedral, atoms
# output - horizonal x id list, vertical x id list
# output - positions corrected, for hor and ver
def find_x_hor_ver(x_oct_found, x_offset_found, atoms):
    x_position_corrected = []
    x_hor_final = []
    x_ver_final = []
    x_hor_position_c = []
    x_ver_position_c = []
    # for corrected position, its a position corrected by periodic boudary
    # x_position_corrected is a list of position in the order of x_oct
    z_position_corrected = []
    for x, off in zip(x_oct_found, x_offset_found):
        x_position_corrected.append(atoms.positions[x] + np.dot(off, atoms.get_cell()))
    for x in range(len(x_oct_found)):
        z_position_corrected.append(x_position_corrected[x][2])
    # z_smallest3 = heapq.nsmallest(3, z_position_corrected)
    # z_avg = max(z_smallest3)
    z_avg = sum(z_position_corrected) / len(z_position_corrected)
    for x in range(len(x_oct_found)):
        if z_avg + 1.5 > x_position_corrected[x][2] > z_avg - 1.5:
            x_hor_final.append(x_oct_found[x])
            x_hor_position_c.append(x_position_corrected[x])
        else:
            x_ver_final.append(x_oct_found[x])
            x_ver_position_c.append((x_position_corrected[x]))
    return x_hor_final, x_ver_final, x_hor_position_c, x_ver_position_c


# test run
if __name__ == "__main__":
    # import numpy as np
    test_atoms = read("../Structs/test/71.vasp", format='vasp')
    # print(len(test_atoms))
    # print(b_finding(test_atoms)[0])
    num_oct = len(b_finding(test_atoms))
    oct_list = []
    for b in b_finding(test_atoms):
        # print("now is b=%s"%b)
        x_oct, x_offset = find_x_oct(b, test_atoms)
        x_hor, x_ver, x_hor_pc, x_ver_pc = find_x_hor_ver(x_oct, x_offset, test_atoms)
        oct_list.append(OctObj(b, x_hor, x_ver, x_hor_pc, x_ver_pc, test_atoms))
        print(b, x_hor, x_ver)
    # b = 15
    # x_oct, x_offset = find_x_oct(b, test_atoms)
    # print(x_oct)
    # print(x_offset)
    # x_hor, x_ver, x_hor_pc, x_ver_pc = find_x_hor_ver(x_oct, x_offset, test_atoms)
    # print(x_ver,x_ver_pc)


    ref_oct = oct_list[0]
    test_oct = oct_list[1]
    from Rotation import rotation_degree

    rot_ang = rotation_degree(ref_oct, test_oct)
    print(rot_ang)

    from Distortion import *

    dis_ang = ang_distortion(ref_oct, test_oct)
    print(dis_ang)
    dis_ver = ver_dist(ref_oct, test_oct)
    print(dis_ver)
    dis_hor = hor_dist(ref_oct, test_oct)
    print(dis_hor)
