import selfies as sf
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from src.objective.guacamol.get_objective import get_objective
from typing import Union
import crossovers as xo

def compound_similarity(smiles, target_smile):
    #print('met.py - compound_similarity input',smiles[0],"\n",smiles)
    fp_target = xo.get_fingerprint(target_smile)
    return [xo.fingerprint_similarity(xo.get_fingerprint(smile),fp_target) for smile in smiles]

def one_over_metric(number):
    if number!=0:
        return 1/number
    return 0

def get_obj(smiles: Union[str,list], objective: str):
    if isinstance(smiles, str): 
        return get_objective([smiles],objective).item()
    return [e[0] for e in get_objective(smiles,objective).tolist()]

def mol_mass_SMILES(smiles: Union[str,list]):
    if isinstance(smiles, str): 
        try:
            return CalcExactMolWt(Chem.MolFromSmiles(smiles))
        except:
            print(f"CalcExactMolWt(Chem.MolFromSmiles(smiles)) err for:\n- SMILES: {smiles}")
            return None
    else:
        return [mol_mass_SMILES(e) for e in smiles]
    
def mol_mass_SELFIES(selfies: Union[str,list]):
    if isinstance(selfies, str): 
        return CalcExactMolWt(Chem.MolFromSmiles(sf.decoder(selfies)))
    return [mol_mass_SELFIES(e) for e in selfies]

def specific_element_count_SMILES(smiles: Union[str,list], element: str):
    if isinstance(smiles, str): 
        return smiles.count(element)
    return [s.count(element) for s in smiles]
