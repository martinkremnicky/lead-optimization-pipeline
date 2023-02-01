import main
import selfies as sf
import random
#import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from molmass import Formula
import mutations as mut
from functools import partial


#just a 'fancy' printing function
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

#print(CalcExactMolWt(Chem.MolFromSmiles(cystene)))
mfl = [partial(mut.mutate,selfies_molecule=sf.encoder(mdma)),
        partial(mut.mutate_insert,selfies_molecule=sf.encoder(mdma),fragment_size=10,random_size=True)]

derivatives = {}

for i in (main.generate_derivatives(100,mfl)):
    print(i)






