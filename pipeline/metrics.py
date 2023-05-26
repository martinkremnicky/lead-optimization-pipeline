import selfies as sf
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from src.objective.guacamol.get_objective import get_objective

def one_over_metric(number):
    if number!=0:
        return 1/number
    return 0

def get_obj(smiles: str, obj: str):
    return get_objective([smiles],obj).item()

def mol_mass_SMILES(smiles: str):
    return CalcExactMolWt(Chem.MolFromSmiles(smiles))
    
def mol_mass_SELFIES(selfies: str):
    return CalcExactMolWt(Chem.MolFromSmiles(sf.decoder(selfies)))

def specific_element_count_SMILES(smiles: str, element: str):
    return smiles.count(element)
