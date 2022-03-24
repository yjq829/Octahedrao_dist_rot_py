# ******************************************************************************************************
# python 3.6
# version 1.0
#
# input part: *.vasp file
# (1) input reference POSCAR
# (2) input test POSCAR
#
# output part:
#
# ******************************************************************************************************

# dependency environment

from ase.io import read


def test_input(test_file):
    file_test = read(test_file, format='vasp')
    return file_test


def get_cell_element(atoms):
    element_list = atoms.get_chemical_symbols()
    element_list = sorted(set(element_list), key=element_list.index)
    # print(element_list)
    return element_list


# test run
if __name__ == "__main__":
    # tot_ele=[]
    # with open("test.csv","w",newline='') as csvfile:
    #     writer=csv.writer(csvfile)
    #     for i in list(range(200)):
    #         test_filename="../Structs/"+str(i+1)+".vasp"
    #         test_atoms=test_input(test_filename)
    #         ele_list=get_cell_element(test_atoms)
    #         tot_ele.extend(ele_list)
    #         writer.writerow([i+1]+ele_list)
    # tot_ele=sorted(set(tot_ele),key=tot_ele.index)
    # print(tot_ele)
    test_filename = "../Structs/1.vasp"
    test_atoms = test_input(test_filename)
    ind_ele = enumerate(test_atoms.get_chemical_symbols())
    b = [i for i, x in ind_ele if x in ["Ge", "Pb"]]
    print(b)
