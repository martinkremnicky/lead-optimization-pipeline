import selfies as sf
import random
#import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from molmass import Formula
import main

alphabet = sf.get_semantic_robust_alphabet() # Gets the alphabet of robust symbols

def mutate(selfies_molecule):
    """
    Generate and swap (SELFIES) molecule's element for another from SELFIES alphabet
    """
    selfies_molecule_list = list(sf.split_selfies(selfies_molecule))
    rnd_index = random.randint(0,len(selfies_molecule_list)-1)
    rnd_sample = random.sample(list(alphabet),1)[0]
    selfies_molecule_list.pop(rnd_index)
    selfies_molecule_list.insert(rnd_index,rnd_sample)
    return main.validate(''.join(selfies_molecule_list))

def mutate_insert(selfies_molecule, fragment_size=3, random_size = True):
    """
    Generate and insert SELFIES fragment into a (SELFIES) molecule 

    :params fragment_size: size of fragment to be generated
        if ``random_size=True``, becomes upper limit of size
        else is exact
    """
    if random_size == True:
        fragment_size = random.randint(1,fragment_size)
    selfies_molecule_list = list(sf.split_selfies(selfies_molecule))
    rnd_index = random.randint(0,len(selfies_molecule_list))
    rnd_sample = random.sample(list(alphabet), fragment_size) 
    rnd_sample = ''.join(rnd_sample)
    selfies_molecule_list.insert(rnd_index,rnd_sample)
    return main.validate(''.join(selfies_molecule_list))