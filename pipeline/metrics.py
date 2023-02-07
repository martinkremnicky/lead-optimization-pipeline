import functions as func
import selfies as sf
import random
#import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from molmass import Formula
import mutations as mut
from functools import partial


def mol_mass_SMILES(smiles_molecule):
    return CalcExactMolWt(Chem.MolFromSmiles(smiles_molecule))
    
def mol_mass_SELFIES(selfies_molecule):
    return CalcExactMolWt(Chem.MolFromSmiles(sf.decoder(selfies_molecule)))

def specific_element_count_SMILES(smiles_molecule: str, element: str):
    return smiles_molecule.count(element)
