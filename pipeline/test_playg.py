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
#import pandas as pd
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

benzene_sf = sf.encoder(benzene)  # [C][=C][C][=C][C][=C][Ring1][=Branch1]


#def test(variable): #TODO a function that prints a var.'s name and value, for debugging
#    print([ i for i, a in locals().items() if a == variable][0],": ",variable)

mutation_function_list = [
    partial(mut.replacement),
    partial(mut.addition,fragment_size=1)
]
metric_function_list = [
    partial(met.mol_mass_SMILES),
    partial(met.specific_element_count_SMILES,element = 'C')
]

n = 200

#pop_df = fn.populate(n,[sf.encoder(mdma),sf.encoder(caffeine)],metric_function_list,mutation_function_list)
pop_df = fn.initialize_pop(n,sf.encoder(caffeine),metric_function_list,mutation_function_list,0,True)


pop_df['Metric 2'] = pop_df['Metric 2'].apply(met.one_over_metric)
#print((fn.classify_fronts(pop_df).tolist()))
print(pop_df)
#for index, row in pop_df.iterrows():
#    print(fn.check_row_dominance(pop_df,index))
dimension_columns = ['Metric 1', 'Metric 2']
#pop_df = fn.pareto_front_classification(pop_df,metrics,minimize=False)
#print(pop_df)
#fn.dominates(pop_df,0)

if False:
    for item in pop_df['SMILES molecule'].tolist():
        print(item)

if False:
    line()
    # these molecule "sets" are the same molecule, but rotated differently
    # (http://www.cheminfo.org/Chemistry/Cheminformatics/Smiles/index.html)
    # trying to figure out how to "derotate" or "match" them/parts of them

    specimen_a_1 = "C1=CC=CC1"
    specimen_a_2 = "C1C=CC=C1"
    print(sf.encoder(specimen_a_1),'\n{}'.format(sf.encoder(specimen_a_2)))

    specimen_b_1 = "C1=CC1"
    specimen_b_2 = "C1C=C1"
    print(sf.encoder(specimen_b_1),'\n{}'.format(sf.encoder(specimen_b_2)))

    specimen_c_1 = "C1=CC=C[S-1]=C1"
    specimen_c_2 = "[S-1]1=CC=CC=C1"
    print(sf.encoder(specimen_c_1),'\n{}'.format(sf.encoder(specimen_c_2)))

    line()
    specimen_d_1 = "C1=[P-1]C=CC=C1"
    specimen_d_2 = "C1=CC=C[P-1]=C1"
    specimen_d_3 = "[P-1]1=CC=CC=C1"
    print(sf.encoder(specimen_d_1),'\n{}\n{}'.format(sf.encoder(specimen_d_2),
                                                        sf.encoder(specimen_d_3)))
    
    # studying SELFIES notation (position of rings and branches)
    line()
    print(sf.encoder("C1=CC=CC=C1[P-1](=[B-1])[N+1]Br"))
    print(sf.encoder(caffeine))
    print(sf.encoder(_12dibromobenzene))
    print(sf.encoder(_13dibromobenzene))
    print(sf.encoder(_14dibromobenzene))

frontier = fn.get_pareto_optimal(pop_df,dimension_columns,minimize=False)
frontier = fn.get_isolation(frontier,dimension_columns)

line(" Molecules on the front: ")
print(frontier)
#print("-/*-**/*--*-*-/-*-*/-*-\n",)

#i_s = (frontier['Isolation score'].tolist())
#sns.scatterplot(data=pop_df, x="Metric 1", y="Metric 2", hue="Isolation score")
#plt.show()

#from collections import Counter
#i_s = [float(i)/sum(i_s) for i in i_s] #normalize
#for i in range(100):
#ct = ((np.random.choice(pop_df.index.tolist(),size=1000,p=i_s))).tolist() #sample 100 from biased roulette wheel
#print(type(ct),ct)
#ct = [str(x) for x in ct] # convert to string for dictionary keys

#print(dict((x,ct.count(x)/10) for x in set(ct))) #count each element
#print([round(x,2) for x in i_s])

line("asssssssss")
if False:
    pop_df2 = fn.populate_list(0,[sf.encoder(mdma),sf.encoder(caffeine)],
                          metric_function_list,mutation_function_list)
#pop_df2.head(100)
#lll = ((pop_df['Origin'].tolist()))
#print(dict((x,lll.count(x)) for x in set(lll)))
if False:
    for a,b in (zip([0,1],[4,5])):
        print(a,b)

frontier_molecules = frontier['SMILES molecule'].tolist()
weights = frontier['Isolation score'].tolist()
frontier_weights_norm = [round(((i)*100)/sum(weights)) for i in weights]
#print(frontier_weights_norm)
import pandas as pd

new_gen = pd.DataFrame()
for molecule, n in zip(frontier_molecules,frontier_weights_norm):
    new_gen = pd.concat([new_gen, fn.initialize_pop(n,sf.encoder(molecule),
                                              metric_function_list,
                                              mutation_function_list,1)])
new_gen.reset_index(drop=True,inplace=True)
new_gen['Metric 2'] = new_gen['Metric 2'].apply(met.one_over_metric)
#print(new_gen)
print(fn.get_pareto_optimal(new_gen,dimension_columns,minimize=False))