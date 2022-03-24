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
# all reference state refer to CsPbCl
# ******************************************************************************************************
# ******************************************************************************************************


# dependency environment

import csv

import IO
import ObjOct
import Rotation
import Distortion

# input reference state
# reference state 1:  CsPbCl

ref_filename_1 = "../Structs/16.vasp"
ref_atoms_1 = IO.test_input(ref_filename_1)
num_ref_oct = len(ObjOct.b_finding(ref_atoms_1))
ref_oct_list = []
for b_ref in ObjOct.b_finding(ref_atoms_1):
    x_oct_ref, x_offset_ref = ObjOct.find_x_oct(b_ref, ref_atoms_1)
    x_hor_ref, x_ver_ref, x_hor_pc_ref, x_ver_pc_ref = ObjOct.find_x_hor_ver(x_oct_ref, x_offset_ref, ref_atoms_1)
    ref_oct_list.append(ObjOct.OctObj(b_ref, x_hor_ref, x_ver_ref, x_hor_pc_ref, x_ver_pc_ref, ref_atoms_1))

default_ref_obj = ref_oct_list[0]

# input test state
# str_input=input("Please input test file number\t")

f = open("../main/result.csv", "w", newline="")
writer = csv.writer(f)
writer.writerow(
    ["No.", "fomula", "Oct1", "", "", "", "Oct2", "", "", "", "Oct3", "", "", "", "Oct4", "", "", "", "Oct5", "", "",
     "", "Oct6", "", "", "", "Oct7", "", "", "", "Oct8", "", "", "", ])
writer.writerow(["", "", "rot_ang", "ang_distortion", "ver_distortion", "hor_distortion", "rot_ang", "ang_distortion",
                 "ver_distortion", "hor_distortion", "rot_ang", "ang_distortion", "ver_distortion", "hor_distortion",
                 "rot_ang", "ang_distortion", "ver_distortion", "hor_distortion", "rot_ang", "ang_distortion",
                 "ver_distortion", "hor_distortion", "rot_ang", "ang_distortion", "ver_distortion", "hor_distortion",
                 "rot_ang", "ang_distortion", "ver_distortion", "hor_distortion", "rot_ang", "ang_distortion",
                 "ver_distortion", "hor_distortion"])
for name in range(441):
    print("Now scanning file # %s" % (name + 1))
    test_filename = "../Structs/" + str(name + 1) + ".vasp"
    test_atoms = IO.test_input(test_filename)
    num_test_oct = len(ObjOct.b_finding(test_atoms))
    test_b_list = ObjOct.b_finding(test_atoms)
    test_oct_list = []
    for b_test in ObjOct.b_finding(test_atoms):
        print("Now working on center B # %s" % b_test)
        x_oct_test, x_offset_test = ObjOct.find_x_oct(b_test, test_atoms)
        x_hor_test, x_ver_test, x_hor_pc_test, x_ver_pc_test = ObjOct.find_x_hor_ver(x_oct_test, x_offset_test,
                                                                              test_atoms)
        test_oct_list.append(ObjOct.OctObj(b_test, x_hor_test, x_ver_test, x_hor_pc_test, x_ver_pc_test, test_atoms))

    # calculate output information
    test_fomula = test_atoms.get_chemical_formula()
    rot_ang_test = []
    dis_ang_test = []
    dis_ver_test = []
    dis_hor_test = []
    result_out = []
    for i in range(num_test_oct):
        rot_ang_test.append(Rotation.rotation_degree(default_ref_obj, test_oct_list[i]))
        dis_ang_test.append(Distortion.ang_distortion(default_ref_obj, test_oct_list[i]))
        dis_ver_test.append((Distortion.ver_dist(default_ref_obj, test_oct_list[i])))
        dis_hor_test.append(Distortion.hor_dist(default_ref_obj, test_oct_list[i]))
    for i in range(num_test_oct):
        result_out.append(rot_ang_test[i])
        result_out.append(dis_ang_test[i])
        result_out.append(dis_ver_test[i])
        result_out.append(dis_hor_test[i])
    writer.writerow([name + 1, test_fomula] + result_out)
f.close()
# output results to screen
# for i in range(num_test_oct):
#     print("******************************************************************************")
#     print("                     The octahedra testing is %s             "%(i+1))
#     print("                     The octahedra fomular is %s             "%test_fomula)
#     print("            The octahedra rotatiion angle is %.3f             "%rot_ang_test[i])
#     print("            The octahedra angular distortion is %.3f             "%dis_ang_test[i])
#     print("            The octahedra vertical distortion is %.3f             "%dis_ver_test[i])
#     print("            The octahedra horizonal distortion is %.3f             "%dis_hor_test[i])
