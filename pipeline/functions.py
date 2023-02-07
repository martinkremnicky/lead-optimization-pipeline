import selfies as sf
import random
#import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from molmass import Formula
import mutations as mut
import numpy as np
import pandas as pd
import copy

alphabet = sf.get_semantic_robust_alphabet() # Gets the alphabet of robust symbols

def validate(selfies_molecule):
    """
    Converts a (generated) SELFIES molecule to SMILES and back 'validate' it
    """
    conversion_smi = sf.decoder(selfies_molecule)
    conversion_sf = sf.encoder(conversion_smi)
    return conversion_smi, conversion_sf

def generate_derivatives(total_size:int,selfies_molecule,mutation_function_list): #add mutation methods parameters
    """
    Creates derivatives of SELFIES molecule, outputs SMILES list by default.
    
    :param total_size: 
    """
    finished_derivates = []
    mfl_len = len(mutation_function_list)
    
    if (total_size % mfl_len != 0):
        pass #raise ... ?
    
    for i in range(total_size):
        finished_derivates.append(mutation_function_list[i%mfl_len](selfies_molecule = selfies_molecule)[0])

    return finished_derivates

def sort_df_by_col_index(df,index: int):
    return df.sort_values(list(df.columns)[index])

def populate(n:int,selfies_molecule,metric_function_list,mutation_function_list):
    """
    Creates a data frame of n derived (SMILES) molecules, with evaluated metrics as additional columns
    """
    column_names = ["SMILES molecule"]
    #for each metric, add "Metric x" as column name
    for i in range(len(metric_function_list)):
        column_names.append("Metric "+str(i+1))
    derivatives = pd.DataFrame(columns=column_names) 
    #while there are fewer than n unique molecules, keep generating
    while derivatives.drop_duplicates().shape[0] < n:
        #for each molecule i
        for i in generate_derivatives(n-derivatives.drop_duplicates().shape[0],selfies_molecule,mutation_function_list):
            #add molecule i to row
            data = [i]
            for j in metric_function_list:
                #add metric evaluations of i
                data.append(j(i))
            #add a row containing the molecule and evaluated metrics to df
            derivatives = pd.concat([derivatives.drop_duplicates(),pd.DataFrame(data=[data],columns=column_names)])  
    derivatives = sort_df_by_col_index(derivatives,1)#sort by 1st metric
    derivatives.index = list(range(derivatives.shape[0]))#reindex
    return derivatives



