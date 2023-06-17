from matplotlib import pyplot as plt
import functions as fn
import selfies as sf
#import random
#import rdkit
#from rdkit import Chem
#from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
#from molmass import Formula
import mutations as mut
from functools import partial
import metrics as met
import pandas as pd
#import copy
import numpy as np
import functions as fn
#import seaborn as sns

#just a 'fancy' printing function to separate outputs
def line(string=''):
    print(f"-------{string}-------")

# smiles
benzene = "c1ccccc1" 
cysteine = "C([C@@H](C(=O)O)N)S"
caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
mdma = "CC(NC)CC1=CC=C(OCO2)C2=C1"
_12dibromobenzene = "Brc1ccccc1Br"
_13dibromobenzene = "Brc1cccc(Br)c1"
_14dibromobenzene = "Brc1ccc(Br)cc1"
_15dibromobenzene = "Brc1cc(Br)ccc1"
_16dibromobenzene = "Brc1c(Br)cccc1"

#study_chem = "CCCCCCCC3CCC(CCCC2CCCC(CCCCCc1cccc(CCCCCC)c1)C2)C3"
#study_chem_sf = sf.encoder(study_chem)
#benzene_sf = sf.encoder(benzene)  # [C][=C][C][=C][C][=C][Ring1][=Branch1]

#print(f"{study_chem}\n{study_chem_sf}")

#def test(variable): #TODO a function that prints a var.'s name and value, for debugging
#    print([ i for i, a in locals().items() if a == variable][0],": ",variable

#print(fn.canonicalize(_12dibromobenzene)+"\n"+fn.canonicalize(_16dibromobenzene))
#print(fn.canonicalize(_13dibromobenzene)+"\n"+fn.canonicalize(_15dibromobenzene))
def generate_derivatives_batch(n:int, selfies_molecule: str, mutation_function_list: list) -> list:



    #"""
    #Creates derivatives of SELFIES molecule, outputs SMILES list by default.
    #"""
    idx = np.arange(n) % len(mutation_function_list)
    print( np.arange(n), len(mutation_function_list),  idx)
    #return [mutation_function_list[i](selfies_molecule)[0] for i in idx]
    

    return ...
three = ['a','b','c']
one = ['a']
generate_derivatives_batch(5,"a",three)