#Author: Phoenix Logan
#Usage: write mmol2 file as sys.argv parameter, must be clean mol2 file
#ex: "python path/to/prgrm only-one-flat.seqn.cln.mol2"

import sys
import pandas as pd
from pandas import DataFrame

def read_mmol2(fh):
    ''' Use pandas to read in mmol2 file and split into a dataframe for
    future parsing'''
    mmol2 = pd.read_table(fh, names = ["col1"])
    mmol2_2 = DataFrame(mmol2)
    return mmol2_2
    
def get_index(identifier, dframe):
    '''break apart mmol2 file by molecule, atom, and bond to parse accordingly
    make sure we are not removing duplicate atoms'''
    #create list of atoms for each molecule after @<TRIPOS>
    atom_list = dframe[dframe["col1"] == "@<TRIPOS>"+identifier].index.tolist()
    if identifier == "MOLECULE":
        atom_list.append(len(dframe.index))
    return atom_list  

def create_chunks(molecule_list, bond_list, atom_list): 
    '''takes the appropriate chunks from mmol2 file (to 
    be parsed and then to be sent to file '''   
    sections = []
    w_sections = []
    for i in range(len(atom_list)):
        range_int = [atom_list[i]+1, bond_list[i]]
        sections.append(range_int) 
    for i in range(len(molecule_list)):
        if i != len(molecule_list)-1:
            range_w = [molecule_list[i], molecule_list[i+1]]
            w_sections.append(range_w)
    return sections, w_sections

def create_dict(section, dframe):
    '''create a dictionary with separated dataframe entries to store each individual molecule'''
    mmol2_dict = {}
    for i in range(len(section)):
        #splice dataframe to get out each molecule
        df_edit = dframe.loc[range(section[i][0], section[i][1])]
        mmol2_dict[i] = df_edit
    return mmol2_dict    

def split(mmol2_dict): 
    ''' splits into columns by white space delimeters using REGEX ''' 
    for i in mmol2_dict.keys():
        j = mmol2_dict.get(i)["col1"].str.split('\s+')
        mmol2_dict[i] = j
    return mmol2_dict
        
def find_flats(mmol2_dict, mmol2_whole): 
    ''' finds the flats, outputs them to new file and extracts them from the mmol2 file '''
    for i in mmol2_dict.keys():
        val = mmol2_dict[i]
        # store z coordinates to determine which molecules have a zero summation"
        z_coords = []
        for k in val:
            z_num = float(k[5])
            z_coords.append(z_num)
        if sum(z_coords) == 0:
            mmol2_whole[i].to_csv("flat_"+str(i+1)+".mol2", index=False, header= False )
            del mmol2_whole[i]
    for value in mmol2_whole.values():
        fh =  open("mmol2_parsed.mol2", "a")
        value.to_csv(fh, mode= "a", header = False, index=False)
        fh.close()
            
def main():
    #run functions
    mmol2_df = read_mmol2(sys.argv[1])  
    molecule_list, bond_list, atom_list = get_index("MOLECULE", mmol2_df), get_index("BOND", mmol2_df), get_index("ATOM", mmol2_df)
    p_sections, w_sections = create_chunks(molecule_list, bond_list, atom_list)
    mmol2_dict_p, mmol2_dict_w = create_dict(p_sections, mmol2_df), create_dict(w_sections, mmol2_df)
    new_mmol2_dict = split(mmol2_dict_p)    
    find_flats(new_mmol2_dict, mmol2_dict_w)
    
if __name__ == "__main__":
    main()
