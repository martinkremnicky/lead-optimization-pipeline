import main
import selfies as sf
import random
#import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from molmass import Formula
import mutations as mut
from functools import partial
import metrics as met
import pandas as pd
import copy
import numpy as np


#just a 'fancy' printing function to separate outputs
def line(string=''):
    print("-------"+string+"-------")

# smiles
benzene = "c1ccccc1" 
cystene = "C([C@@H](C(=O)O)N)S"
caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
mdma = "CC(NC)CC1=CC=C(OCO2)C2=C1"

benzene_sf = sf.encoder(benzene)  # [C][=C][C][=C][C][=C][Ring1][=Branch1]


#def test(variable): #TODO
#    print([ i for i, a in locals().items() if a == variable][0],": ",variable)

mutation_function_list = [
    partial(mut.mutate),
    partial(mut.mutate_insert)
]
metric_function_list = [
    partial(met.mol_mass_SMILES),
    partial(met.specific_element_count_SMILES,element = 'C')
]

print(main.populate(100,benzene_sf,metric_function_list,mutation_function_list))
