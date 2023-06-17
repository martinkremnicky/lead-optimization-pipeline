import selfies as sf
import numpy as np
import pandas as pd
from rdkit import Chem
from typing import Tuple, Union, List
import random
from constants import *
import functions as fn
import mutations as mut
import metrics as met

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

from rdkit import DataStructs

import time


def random_individual_crossover(selfies_list: List[str],n:int, crossover_type = 0) -> List[str]:
    """
    Returns `n` offspring of two randomly selected parents.
    """
    new_pool = []
    for _ in range(n):
        j = random.randint(0,len(selfies_list)-1)
        k = random.randint(0,len(selfies_list)-1)
        if crossover_type == 0:
            new_pool.append(get_path_middle(selfies_list[j], selfies_list[k]))
        elif crossover_type == 1:
            new_pool.append(get_median_molecule(selfies_list[j], selfies_list[k]))
        else:
            print("Incorrect crossover type selected")
    return new_pool

def get_median_molecule(selfies1: str, selfies2: str) -> str:
    #print(f"get_median_molecule({selfies1}, {selfies2})")
    chem_path_smiles = [str(sf.decoder(chem)) for chem in chemical_path(selfies1, selfies2)+chemical_path(selfies2, selfies1)]

    smiles1 = str(sf.decoder(selfies1))
    smiles2 = str(sf.decoder(selfies2))

    joint_sim_list = []
    for chem in chem_path_smiles: 
        joint_sim_list.append(joint_sim([smiles1,smiles2],chem))
    max_js =  max(joint_sim_list) #biggest JS
    return sf.encoder(chem_path_smiles[joint_sim_list.index(max_js)])



def joint_sim(M: List[str],m: str) -> float:
    sim_sum = 0
    sim_list = []
    for m_i in M:
        sim_sum += fingerprint_similarity(m_i,m)
        sim_list.append(fingerprint_similarity(m_i,m))
    sim_sum = sim_sum/len(M)
    return sim_sum - (max(sim_list) - min(sim_list))

def fingerprint_similarity(smiles1:str, smiles2:str) -> float:
    #print(f"fingerprint_similarity({smiles1},{smiles2})")
    mols = [Chem.MolFromSmiles(smiles1), Chem.MolFromSmiles(smiles2)]
    fpgen = AllChem.GetRDKitFPGenerator()
    fps = [fpgen.GetFingerprint(x) for x in mols]
    return DataStructs.TanimotoSimilarity(fps[0], fps[1])




def get_path_middle(selfies1:str, selfies2:str) -> str: 
    """returns SELFIES """
    chem_path = chemical_path(selfies1, selfies2)
    # 0 1 2 3 4 5 6 -> 3 = 6/2
    # 0 1 2 3 4 5  -> 2,3 = floor(5/2), round(5/2) #actually ill just round so I pick one offspring
    return chem_path[round(len(chem_path)/2)]

def chemical_path(selfies1: str, selfies2: str) -> List[str]:
    """returns list of SELFIES"""
    sf_lst1 = list(sf.split_selfies(selfies1))
    sf_lst2 = list(sf.split_selfies(selfies2))

    sf_lst1_len = len(sf_lst1)
    sf_lst2_len = len(sf_lst2)
    
    if sf_lst1_len>sf_lst2_len:
        sf_lst2.extend(['[nop]']*(sf_lst1_len-sf_lst2_len))
    elif sf_lst2_len>sf_lst1_len:
        sf_lst1.extend(['[nop]']*(sf_lst2_len-sf_lst1_len))

    chem_path = [(selfies1)]
    for sym_i in range(sf_lst1_len):
        sf_lst1[sym_i] = sf_lst2[sym_i]
        chem_path.append((''.join(sf_lst1)))

    return chem_path

