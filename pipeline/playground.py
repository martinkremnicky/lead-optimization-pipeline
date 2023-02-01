import main
import selfies as sf
import random
#import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from molmass import Formula
import mutations as mt
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

def func(a,b):
    return random.randint(a,b)

#a_func = func(1,5)
func_list = [partial(func,a=1,b=5), partial(func,a=1,b=5), partial(func,a=1,b=5)]

print("output:")
for i in range(4):
    for j in range(len(func_list)):
        print(func_list[j]())

#print(CalcExactMolWt(Chem.MolFromSmiles(cystene)))
mfl = [partial(mt.mutate,selfies_molecule=benzene_sf),
        partial(mt.mutate_insert,selfies_molecule=benzene_sf,fragment_size=10,random_size=True)]

derivatives = {}

for i in (main.generate_derivatives(100,mfl)):
    print(i)






