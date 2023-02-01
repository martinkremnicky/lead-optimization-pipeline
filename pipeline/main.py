import selfies as sf
import random
#import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from molmass import Formula
import mutations as mt

alphabet = sf.get_semantic_robust_alphabet() # Gets the alphabet of robust symbols

def validate(selfies_molecule):
    """
    Converts a (generated) SELFIES molecule to SMILES and back 'validate' it
    """
    conversion_smi = sf.decoder(selfies_molecule)
    conversion_sf = sf.encoder(conversion_smi)
    return conversion_smi, conversion_sf

def generate_derivatives(total_size,mutation_function_list): #add mutation methods parameters
    """
    Creates derivatives of SELFIES molecule, outputs SMILES list by default.
    
    :param total_size: 
    """
    finished_derivates = []#[sf.decoder(selfies_molecule)] #add the original in case of population-based algorithms
    for i in range(int(total_size/(len(mutation_function_list)-1))):
        print('i:',i)
        for j in range(len(mutation_function_list)):
            finished_derivates.append(mutation_function_list[j]()[0])
    return finished_derivates

def mol_mass_SMILES(smiles_molecule):
    return CalcExactMolWt(Chem.MolFromSmiles(smiles_molecule))
    
def mol_mass_SELFIES(selfies_molecule):
    selfies_molecule = sf.decoder(validate(selfies_molecule))
    return CalcExactMolWt(Chem.MolFromSmiles(selfies_molecule))


