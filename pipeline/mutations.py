import selfies as sf
import random
import functions as fn
import numpy as np
import copy
from constants import *


def replacement(selfie, n=1, random_size = False, ring_aware = False): #TODO: add 'n' swaps?
    """
    Generate and swap (SELFIES) molecule's element for another from SELFIES alphabet
    """
    min_size = 2 #hyperparam to prevent errors (empty SMILES) #1
    try:
        selfies_molecule_list = list(sf.split_selfies(selfie))
    except:
        print(f"split_selfies err for input: (check if SELFIES or SMILES)\n{selfie}")
    selfies_size = len(selfies_molecule_list)
    if selfies_size>min_size:
        rnd_sample = random.sample(list(ALPHABET),1)[0]
    else:
        rnd_sample = random.sample(list(TRANSLATABLE_ALPHABET),1)[0]
    if random_size:
        n = random.randint(1,n)
    for _ in range(n):
        rnd_index = random.randint(0,selfies_size-1)
        if ring_aware:
            try:
                if selfies_molecule_list[rnd_index-1] in ['[Ring1]','[Ring2]']:
                    rnd_sample = random.sample(['[Ring1]','[Ring2]','[Branch1]','[=Branch1]'],1)[0]
            except:
                pass
        selfies_molecule_list.pop(rnd_index)
        selfies_molecule_list.insert(rnd_index,rnd_sample)
    try:
        return fn.validate(''.join(selfies_molecule_list))
    except:
        replacement(selfie, n, random_size, ring_aware)
    
def addition(selfies_molecule, fragment_size=1, random_size = False, rings = True):
    """
    Generate and insert SELFIES fragment into a (SELFIES) molecule 

    ### Parameters 
    fragment_size: size of fragment to be generated
        - if ``random_size=True``, becomes upper size limit of fragment generated,
        else (if `=False`) becomes exact fragment size
    """
    if random_size:
        fragment_size = random.randint(1,fragment_size)
    selfies_molecule_list = list(sf.split_selfies(selfies_molecule))
    rnd_index = random.randint(0,len(selfies_molecule_list))
    pool = list(ALPHABET)
    if rings:
        pool += [sf.encoder(i) for i in ['C1CC1','C1CCC1','C1CCCC1','C1CCCCC1','c1ccccc1']]

    rnd_sample = random.sample(pool, fragment_size) 
    rnd_sample = ''.join(rnd_sample)
    selfies_molecule_list.insert(rnd_index,rnd_sample)
    try:
        return fn.validate(''.join(selfies_molecule_list))
    except:
        addition(selfies_molecule,fragment_size,random_size,rings)

def deletion(selfies_molecule, n=1, random_size = False): 
    selfies_molecule_list = list(sf.split_selfies(selfies_molecule))
    if len(selfies_molecule_list) == 1:
        return fn.validate(''.join(selfies_molecule_list))
    if random_size:
        n = random.randint(1,n)
    if n >= len(selfies_molecule_list):
        #raise ValueError(f"'n' should not be greater than or equal to the length of the molecule (n={n} vs len(mol)={len(selfies_molecule_list)} for {selfies_molecule})")
        n-=len(selfies_molecule_list)-1
    idx = random.sample(range(len(selfies_molecule_list)), n)
    
    # sort indices in descending order so deleting doesn't mess up subsequent indices
    idx.sort(reverse=True)

    for i in idx:
        del selfies_molecule_list[i]
    try:
        return fn.validate(''.join(selfies_molecule_list))
    except:
        deletion(selfies_molecule, n, random_size)

def none(selfies_molecule):
    return fn.validate(selfies_molecule)