#!/usr/bin/python


#usage: script.py  xtal_lig  recp



from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os
from pandas import DataFrame
import glob

def build_lig_file_list():
    pdbqt_file_list = glob.glob("*.pdbqt")
    #pdbqt_file_list.extend(glob.glob("*0001/ZINC*.pdbqt"))
    mol2_file_list = glob.glob("*.mol2")
    #mol2_file_list.extend(glob.glob("*0001/ZINC*.mol2"))
    return pdbqt_file_list, mol2_file_list

def extract_atom_types( ligand_file_list ):
    ''' extract the unique atom types for AD4 gpf file'''
    # list all of the pdbqt files to find unique atom types
    file_lines = []
    # open each file
    for i in ligand_file_list:
        j = open(i,"r")
        for line in j:
            # split each line by its space
            k = line.split(" ")
            sorted_k = []
            for l in k:
                # remove empty list elements
                if len(l) > 0:
                    sorted_k.append(l)
            file_lines.append(sorted_k)
    # turn list of file lines into a dataframe to make extraction easier
    frame = DataFrame(file_lines)
    # drop the duplicates
    type_lst = frame[12].drop_duplicates()
    new_types = []
    # remove atom types if they are NoneType or \n
    for i in type_lst:
        try:
            if i[-1] == "\n":
                j = i[:len(i)-1]
                new_types.append(j)
            else:
                new_types.append(i)
        except TypeError:
            pass
    return new_types

def dock_box( molfile ):
    # compute center coordinates (vector) and box edge lengths (vector)
    # read in the crystal ligand and extract the centroid 
    # for finding the correct box dimensions
    m = Chem.MolFromMol2File( molfile, sanitize = False )

    # center coords
    conf1 = m.GetConformer()
    centroid = Chem.rdMolTransforms.ComputeCentroid(conf1)
    x_center = centroid.x
    y_center = centroid.y
    z_center = centroid.z

    # get molecules coords
    atom_coord_x = []
    atom_coord_y = []
    atom_coord_z = []
    for atom in m.GetAtoms():
        atom_index = atom.GetIdx()
        pos = conf1.GetAtomPosition(atom_index)
        atom_coord_x.append(pos.x)
        atom_coord_y.append(pos.y)
        atom_coord_z.append(pos.z)
    # calc box dimensions basd on min/max of xyz and pad with 5 angstroms
    padding = 4
    x_length = (max(atom_coord_x) + padding) - (min(atom_coord_x) - padding)
    y_length = (max(atom_coord_y) + padding) - (min(atom_coord_y) - padding)
    z_length = (max(atom_coord_z) + padding) - (min(atom_coord_z) - padding)

    box_center = (x_center, y_center, z_center)
    box_dim = (x_length, y_length, z_length)
    return box_center, box_dim

def add_uniques( gpf, atm_types_lst, receptor_type_lst, box_center, box_dim, recp ):
    # this function will add the unique elements to the gpf file
    # to create new grid maps
    os.chdir("shared")
    lig_str_types = ""
    recp_str_types = ""
    map_types = ""
    for i in atm_types_lst:
        lig_str_types += i+" "
        j = "map receptor."+i+".map                   # atom-specific affinity map\n"
        map_types += j
    for i in receptor_type_lst:
        recp_str_types += i+" "
    print "unique"+lig_str_types
    print "unique"+recp_str_types
    # modify the gpf file
    # 40.0 * 0.375 = 15.0 Angstrom box edge lengths...
    x_len, y_len, z_len = [ str(int(round(float(edge_length)/0.375))) for edge_length in box_dim ]
    # must consider what smina has been doing based on crystal ligand dimensions
    with open(gpf, "w") as fh:
        fh.write(
    "npts  "+x_len+" "+y_len+" " +z_len+" # num.grid points in xyz\n"
    "gridfld receptor.maps.fld            # grid_data_file\n"
    "spacing 0.375                        # spacing(A)\n"
    "receptor_types "+recp_str_types+"    # receptor atom types\n"
    "ligand_types "  +lig_str_types+"     # ligand atom types\n"
    "receptor "      + recp +  "          # macromolecule\n"
    "gridcenter "+' '.join([str(x) for x in box_center])+" # xyz-coordinates or auto\n"
    "smooth 0.5                           # store minimum energy w/in rad(A)\n"
    +map_types+
    "elecmap receptor.e.map               # electrostatic potential map\n"
    "dsolvmap receptor.d.map              # desolvation potential map\n"
    "dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, cons\n"
    )
    
    os.system("autogrid4 -p "+gpf+" -l recep.glg")



def AutoDock4(mol_lig, pdbqt_lig, u_lst):
    ''' This function will write out Autodock4 docking parameter file (dpf) '''
    try:
        num_tor_file = open(pdbqt_lig, "r")
        num_tors_lst = num_tor_file.readlines()[0].split(" ")
        num_tors = str(num_tors_lst[2])
    except IndexError:
        print "BAD PDBQT"

    # for building ligand atom-type interaction maps for the receptor
    str_types = ""
    str_maps = ""
    for i in u_lst:
       str_types += i+" "
       j = "map receptor."+i+".map                   # atom-specific affinity map \n"
       str_maps += j
    # try except statements to keep program running through 'lemon' molecules
    try:
        m = Chem.MolFromMol2File(mol_lig, sanitize = False)
        conf1 = m.GetConformer()
        centroid = Chem.rdMolTransforms.ComputeCentroid(conf1)
        x_center = float(str(centroid.x)[:4])
        y_center = float(str(centroid.y)[:4])
        z_center = float(str(centroid.z)[:4])

        config_file = open(pdbqt_lig[:len(pdbqt_lig)-6]+"_config.dpf","w")
        config_file.write(
        "autodock_parameter_version 4.2       # used by autodock to validate parameter set \n"

        "outlev 1                             # diagnostic output level \n"
        "intelec                              # calculate internal electrostatics \n"
        "seed pid time                       # seeds for random generator \n"
        "ligand_types "+str_types+"             # atoms types in ligand \n"
        "fld receptor.maps.fld                # grid_data_file \n"
        
        +str_maps+

        "elecmap receptor.e.map               # electrostatics map \n"
        "desolvmap receptor.d.map             # desolvation map \n"
        "move "+pdbqt_lig+" # small molecule \n"
        "about "+str(x_center)+" "+str(y_center)+" "+str(z_center)+" # small molecule center \n"
        "tran0 random                         # initial coordinates/A or random \n"
        "quaternion0 random                   # initial orientation \n"
        "dihe0 random                         # initial dihedrals (relative) or random \n"
        "torsdof "+num_tors+"             # torsional degrees of freedom \n"
        "rmstol 2.0                           # cluster_tolerance/A \n"
        "extnrg 1000.0                        # external grid energy \n"
        "e0max 0.0 10000                      # max initial energy; max number of retries \n"
        "ga_pop_size 150                      # number of individuals in population \n"
        "ga_num_evals 2500000                 # maximum number of energy evaluations \n"
        "ga_num_generations 27000             # maximum number of generations \n"
        "ga_elitism 1                         # number of top individuals to survive to next ge \n"
        "ga_mutation_rate 0.02                # rate of gene mutation \n"
        "ga_crossover_rate 0.8                # rate of crossover \n"
        "ga_window_size 10                    # \n"
        "ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution \n"
        "ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution \n"
        "set_ga                               # set the above parameters for GA or LGA \n" 
        "sw_max_its 300                       # iterations of Solis & Wets local search \n"
        "sw_max_succ 4                        # consecutive successes before changing rho \n"
        "sw_max_fail 4                        # consecutive failures before changing rho \n"
        "sw_rho 1.0                           # size of local search space to sample \n"
        "sw_lb_rho 0.01                       # lower bound on rho \n"
        "ls_search_freq 0.06                  # probability of performing local search on indiv \n"
        "set_psw1                             # set the above pseudo-Solis & Wets parameters \n"
        "unbound_model bound                  # state of unbound ligand \n"
        "ga_run 10                            # do this many hybrid GA-LS runs \n"
        "analysis                             # perform a ranked cluster analysis"
         )
        config_file.close()
    except (AttributeError, UnboundLocalError) as e:
        print "molecule is lemon"
        pass



def config_sort(total_uniques):
    dir_list = glob.glob("dec*")
    dir_list.extend(glob.glob("act*"))
    for i in dir_list:
        if os.path.isdir(i):
            os.chdir(i)
            pdbqt_lst, mol2_lst = build_lig_file_list()
            uniq_atoms = extract_atom_types( pdbqt_lst )
            for j,k in zip(pdbqt_lst, mol2_lst):
                uniq_atoms = extract_atom_types( pdbqt_lst )
                AutoDock4(k,j,total_uniques)
        os.chdir("..")





def absolute_uniques(gpf_name, box_center, box_dim, recp):
    dir_list = glob.glob("dec*")
    dir_list.extend(glob.glob("act*"))
    total_uniques = []
    recp_uniques = []
    for i in dir_list:
        if os.path.isdir(i):
            os.chdir(i)
            pdbqt_lst, mol2_lst = build_lig_file_list()
            for j,k in zip(pdbqt_lst, mol2_lst):
                uniq_atoms = extract_atom_types( pdbqt_lst )
                for i in uniq_atoms:
                    total_uniques.append(i)
        os.chdir("..")
    uniq_recp = extract_atom_types([recp])
    frame = DataFrame(total_uniques)
    type_lst = frame[0].drop_duplicates()
    new_types = []
    for i in type_lst:
        try:
            if i[-1] == "\n":
                j = i[:len(i)-1]
                new_types.append(j)
            else:
                new_types.append(i)
        except TypeError:
            pass
    print new_types
    return new_types

def main():
    xtal_lig = sys.argv[1]
    recp     = sys.argv[2]

    box_center, box_dim = dock_box(xtal_lig)
    gpf_name = "test.gpf"
    uniques = absolute_uniques(gpf_name, box_center, box_dim, recp)
    config_sort(uniques)
    receptor_type_lst = extract_atom_types([recp])
    add_uniques(gpf_name, uniques, receptor_type_lst, box_center, box_dim, recp)
    

if __name__ == '__main__':
    main()

