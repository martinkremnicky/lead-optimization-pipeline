import selfies as sf
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt

def one_over_metric(number):
    if number!=0:
        return 1/number
    return 0

def mol_mass_SMILES(smiles_molecule):
    return CalcExactMolWt(Chem.MolFromSmiles(smiles_molecule))
    
def mol_mass_SELFIES(selfies_molecule):
    return CalcExactMolWt(Chem.MolFromSmiles(sf.decoder(selfies_molecule)))

def specific_element_count_SMILES(smiles_molecule: str, element: str):
    return smiles_molecule.count(element)
