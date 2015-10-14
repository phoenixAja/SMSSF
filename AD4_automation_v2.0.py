__author__ = 'plogan'


#!/usr/bin/python

'''
Description: script automates preparation of AutoDock4 docking
runs. It first prepares an AutoGrid4 input file and runs AutoGrid
to generate interaction maps for the target. It also makes a list
of all ligand atoms types and creates independent AutoDock4
configuration files for each ligand-target docking.

usage: script.py  targ_name   path_to_xtal_lig   path_to_recp

./AD4_automation_v1.11.py adrb1 ./shared-adrb1/crystal_ligand.mol2 ./shared-adrb1/recp_adrb1.pdbqt
'''

from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os
from pandas import DataFrame
import glob

class Ligand(object):
    def __init__(self, xtal_name):
        self.xtal_name = xtal_name
        self.mol2 = Chem.MolFromMol2File(self, sanitize = False )
        self.conf1 = self.mol2.GetConformer
        self.centroid = Chem.rdMolTransforms.ComputeCentroid(self.conf1)
        self.x_centroid = centroid.x
        self.y_centroid = centroid.y
        self.z_centroid = centroid.z
        self.atom_coord_x = []
        self.atom_coord_y = []
        self.atom_coord_z = []
        for atom in self.mol2.GetAtoms():
            atom_index = atom.GetIdx()
            pos = conf1.GetAtomPosition(atom_index)
            self.atom_coord_x.append(pos.x)
            self.atom_coord_y.append(pos.y)
            self.atom_coord_z.append(pos.z)
        '''compute box center coordinates (vector) and box edge
       lengths (vector). Use crystal ligand's centroid and
       pad min max coords in each dimension for box dimensions'''
        padding = 4
        self.x_length = ( max(self.atom_coord_x) + padding) - (min(self.atom_coord_x) - padding )
        self.y_length = ( max(self.atom_coord_y) + padding) - (min(self.atom_coord_y) - padding )
        self.z_length = ( max(self.atom_coord_z) + padding) - (min(self.atom_coord_z) - padding )
        self.box_dim = (self.x_length, self.y_length, self.z_length)
        self.box_center = (self.x_centroid, self.y_centroid, self.z_centroid)

targ_name = sys.argv[1]
xtal_lig  = Ligand(sys.argv[2])
recp      = sys.argv[3]
gpf_name  = "recp_"+targ_name+".gpf"

def build_lig_file_list():
    '''get mol2 and pdbqt lists of the ligand files'''
    pdbqt_file_list = glob.glob("*.pdbqt")
    mol2_file_list = glob.glob("*.mol2")
    return pdbqt_file_list, mol2_file_list


def extract_atom_types( ligand_file_list ):
    ''' extract the unique atom types for AD4 gpf file'''
    file_lines = []
    for i in ligand_file_list:
        j = open(i,"r")
        for line in j:
            k = line.split(" ")
            sorted_k = []
            for l in k:
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


def build_grids( gpf, atm_types_lst, receptor_type_lst, box_center, box_dim, recp ):
    '''this function will add the unique elements to the gpf file
       and then run AutoGrid to get the grid maps'''
    recp_file = os.path.basename(recp)
    os.chdir("shared-"+targ_name)
    lig_str_types = ""
    recp_str_types = ""
    map_types = ""
    for i in atm_types_lst:
        lig_str_types += i+" "
        j = "map recp_"+targ_name+"."+i+".map                   # atom-specific affinity map\n"
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
    "gridfld recp_"+targ_name+".maps.fld            # grid_data_file\n"
    "spacing 0.375                        # spacing(A)\n"
    "receptor_types "+recp_str_types+"    # receptor atom types\n"
    "ligand_types "  +lig_str_types+"     # ligand atom types\n"
    "receptor "      + recp_file +  "          # macromolecule\n"
    "gridcenter "+' '.join([str(x) for x in box_center])+" # xyz-coordinates or auto\n"
    "smooth 0.5                           # store minimum energy w/in rad(A)\n"
    +map_types+
    "elecmap recp_"+targ_name+".e.map               # electrostatic potential map\n"
    "dsolvmap recp_"+targ_name+".d.map              # desolvation potential map\n"
    "dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, cons\n"
    )

    os.system("autogrid4 -p "+gpf+" -l recp_"+targ_name+".glg")


def get_tors_num( pdbqt_lig ):
    '''read through pdbqt and find number of torsions'''
    with open( pdbqt_lig, 'r') as fh:
        data = fh.readlines()
    for line in data:
        if line.split()[2] == 'active':
            return str(line.split()[1])


def write_AD4_config(mol_lig, pdbqt_lig, u_lst):
    ''' This function will write out Autodock4 docking parameter file (dpf) '''
    try:
        num_tors = get_tors_num( pdbqt_lig )
    except IndexError:
        print "BAD PDBQT"

    # for building ligand atom-type interaction maps for the receptor
    str_types = ""
    str_maps = ""
    for i in u_lst:
       str_types += i+" "
       j = "map recp_"+targ_name+"."+i+".map                   # atom-specific affinity map \n"
       str_maps += j
    try:
        # get centroid for ligand you want to dock, should be present as mol2 file
        Lig = Ligand(mol_lig)
        x_center = float(str(Lig.x_centroid)[:4])
        y_center = float(str(Lig.y_centroid)[:4])
        z_center = float(str(Lig.z_centroid)[:4])


        config_file = open(pdbqt_lig[:len(pdbqt_lig)-6]+"_config_"+targ_name+".dpf","w")
        config_file.write(
        "autodock_parameter_version 4.2       # used by autodock to validate parameter set \n"

        "outlev 1                             # diagnostic output level \n"
        "intelec                              # calculate internal electrostatics \n"
        "seed pid time                       # seeds for random generator \n"
        "ligand_types "+str_types+"             # atoms types in ligand \n"
        "fld recp_"+targ_name+".maps.fld                # grid_data_file \n"

        +str_maps+

        "elecmap recp_"+targ_name+".e.map               # electrostatics map \n"
        "desolvmap recp_"+targ_name+".d.map             # desolvation map \n"
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


def batch_write_configs(total_uniques):
    '''loop through lig dirs and generate AD4 configuration files'''
    dir_list = glob.glob( targ_name+"-*" )
    for ligdir in dir_list:
        if os.path.isdir(ligdir):
            os.chdir(ligdir)
            pdbqt_lst, mol2_lst = build_lig_file_list()
            for i in range(0,len(pdbqt_lst)):
                write_AD4_config( mol2_lst[i], pdbqt_lst[i], total_uniques )
            os.chdir("..")


def absolute_uniques():
    '''get list of unique atoms over the set of ligands'''
    dir_list = glob.glob( targ_name+"-*" )
    total_uniques = []
    for ligdir in dir_list:
        if os.path.isdir(ligdir):
            os.chdir(ligdir)
            pdbqt_lst, mol2_lst = build_lig_file_list()
            uniq_atoms = extract_atom_types( pdbqt_lst )
            for atm in uniq_atoms:
                total_uniques.append(atm)
            os.chdir("..")
    new_types = list(set(total_uniques))
    return new_types


def main():
    uniques = absolute_uniques()
    batch_write_configs(uniques)
    receptor_type_lst = extract_atom_types([recp])
    box_center, box_dim = xtal_lig.box_center, xtal_lig.box_dim
    build_grids(gpf_name, uniques, receptor_type_lst, box_center, box_dim, recp)


if __name__ == '__main__':
    main()

