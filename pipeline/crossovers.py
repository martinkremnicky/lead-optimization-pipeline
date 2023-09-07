import selfies as sf
import numpy as np
import pandas as pd
from rdkit import Chem
from typing import Tuple, Union, List
import random
from constants import *
import functions as fn
import mutations as mut
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import DataStructs
import time
from functools import lru_cache


def joint_sim(M: List[str],m: str) -> float:
    sim_list = np.array([fingerprint_similarity(m_i,m) for m_i in M])
    sim_sum = np.mean(sim_list)
    return sim_sum - (sim_list.max() - sim_list.min())

def fingerprint_similarity(fp1, fp2) -> float:
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def get_fingerprint(smile):
    mol = Chem.MolFromSmiles(smile)
    fpgen = AllChem.GetRDKitFPGenerator()
    return fpgen.GetFingerprint(mol)

def get_median_molecule(selfie1: str, selfie2: str) -> str:
    chem_path_selfies = bidirectional_chem_path(selfie1, selfie2)
    chem_path_smiles = [str(sf.decoder(chem)) for chem in chem_path_selfies]

    smiles1 = str(sf.decoder(selfie1))
    smiles2 = str(sf.decoder(selfie2))
    
    fp1 = get_fingerprint(smiles1)
    fp2 = get_fingerprint(smiles2)

    joint_sim_list = []
    for chem in chem_path_smiles: 
        fp3 = get_fingerprint(chem)
        joint_sim_list.append(joint_sim([fp1,fp2], fp3))
    max_js =  max(joint_sim_list) #biggest JS
    return chem_path_selfies[joint_sim_list.index(max_js)]

def random_individual_crossover(selfies: List[str],n:int, crossover_type = 0) -> List[str]:
    """
    Returns `n` offspring of two randomly selected parents.
    """
    new_pool = []

    selfies = [sf.encoder(fn.canonicalize_selfie(selfie)) for selfie in selfies]
    
    for _ in range(n):
        j = random.randint(0,len(selfies)-1)
        k = random.randint(0,len(selfies)-1)
        if selfies[j] == selfies[k]: #one more attempt at picking two different  individuals
            k = random.randint(0,len(selfies)-1) 



        chem_path = chemical_path(selfies[j], selfies[k])


        
        if crossover_type == 0:
            new_pool.append(get_path_random(chem_path))
        elif crossover_type == 1:
            new_pool.append(get_median_molecule(selfies[j], selfies[k]))
        else:
            print("Incorrect crossover type selected")
    return new_pool


def get_path_random(path):
    i = random.randint(0,len(path)-1)
    return path[i]


def get_path_middle(selfie1:str, selfie2:str) -> str: 
    """returns SELFIES """
    chem_path = chemical_path(selfie1, selfie2)
    return chem_path[round(len(chem_path)/2)]

def chemical_path(selfie1: str, selfie2: str) -> List[str]:
    """returns list of SELFIES"""
    starting_selfie_chars = list(sf.split_selfies(selfie1))
    target_selfie_chars = list(sf.split_selfies(selfie2))
    
    if len(starting_selfie_chars)>len(target_selfie_chars):
        target_selfie_chars.extend(['[nop]']*(len(starting_selfie_chars)-len(target_selfie_chars)))
    elif len(starting_selfie_chars)<len(target_selfie_chars):
        starting_selfie_chars.extend(['[nop]']*(len(target_selfie_chars)-len(starting_selfie_chars)))

    indices_diff = [i for i in range(len(starting_selfie_chars)) if starting_selfie_chars[i] != target_selfie_chars[i]]
    random.shuffle(indices_diff)

    chem_path = [selfie1]
    for sym_i in indices_diff:
        starting_selfie_chars[sym_i] = target_selfie_chars[sym_i]
        chem_path.append((''.join(starting_selfie_chars)))

    return chem_path

def canonical_chem_path(selfie1: str, selfie2: str) -> List[str]:
    """returns list of SELFIES"""
    sf1 = sf.encoder(fn.canonicalize_selfie(selfie1))
    sf2 = sf.encoder(fn.canonicalize_selfie(selfie2))
    return chemical_path(sf1,sf2)

def bidirectional_chem_path(selfie1: str, selfie2: str) -> List[str]:
    return chemical_path(selfie1,selfie2) + chemical_path(selfie2,selfie1)

def randomized_chem_path(selfie1: str, selfie2: str) -> List[str]:
    rs1 = fn.randomize_selfie(selfie1)
    rs2 = fn.randomize_selfie(selfie2)
    return chemical_path(rs1,rs2)

def all_chem_paths(selfie1: str, selfie2: str, n_random = 1) -> List[str]:
    bidir = bidirectional_chem_path(selfie1, selfie2)
    #bidir = []
    canonical = canonical_chem_path(selfie1, selfie2) + canonical_chem_path(selfie2, selfie1)
    rand = []
    if n_random>0:
        for _ in range(n_random):
            rand.extend(randomized_chem_path(selfie1, selfie2))
            rand.extend(randomized_chem_path(selfie2, selfie1))
    #print(len(list(set(bidir + canonical + rand)))/len(list(set(bidir + canonical))))
    return list(set(bidir + canonical + rand))