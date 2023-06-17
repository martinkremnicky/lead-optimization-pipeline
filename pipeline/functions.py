import selfies as sf
import numpy as np
import pandas as pd
from rdkit import Chem
from typing import Tuple, Union, List
import random
from constants import *
import crossovers as xo
from multiprocessing import Pool
import time



def canonicalize_smiles(smiles: str) -> str:
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)
    except:
        print(f"Failed to 'canonicalize_smiles' for {smiles}")
        pass

def canonicalize_selfies(selfies: str) -> str:
    try: 
        return Chem.MolToSmiles(Chem.MolFromSmiles(sf.decoder(selfies)), canonical=True)
    except:
        print(f"Failed to 'canonicalize_selfies' for {selfies}")
        pass

def validate(selfies_molecule:str)->Tuple[str, str]:
    """
    Converts a (generated) SELFIES molecule to SMILES and back to 'validate' it
    ### Output:
    tuple of (SMILES,SELFIES) notations of the same molecule
    """
    conversion_smi = sf.decoder(selfies_molecule)
    try:
        conversion_sf = sf.encoder(conversion_smi)
        return conversion_smi, conversion_sf
    except:
        print(f"functions.py - conversion_sf of smiles ({conversion_smi}) resulted in err (ORIGINAL SELFIES: {selfies_molecule})")

def sort_df_by_col_index(df: pd.DataFrame, index: int) -> pd.DataFrame:
    return df.sort_values(df.columns[index]).reset_index(drop=True)

def column_to_list(df: pd.DataFrame, col_name:str):
    if col_name not in df.columns:
        raise ValueError(f"The DataFrame does not contain a column named '{col_name}'.")
    return df[col_name].tolist()

def get_smiles_list(df: pd.DataFrame, verify=True) -> list:
    accepted_column_names = ["SMILES molecule","smiles", "SMILES"]
    for name in accepted_column_names:
        try:
            if verify:
                return column_to_list(process_df(df,name),name)
            else:
                return column_to_list(df,name)
        except:
            pass
    print(f"Could not find a column with an accepted name. Accepted names:\n{accepted_column_names}") #TODO make `raise`?

def generate_derivatives(n:int, selfies_molecule: str, mutation_function_list: list) -> list:
    """
    Creates `n` derivatives of SELFIES molecule, outputs SMILES list by default.
    """
    idx = np.arange(n) % len(mutation_function_list) #e.g. [0 1 2 3 4 5...] % 3 -> [0 1 2 0 1 2 ...] 
    return [mutation_function_list[i](selfies_molecule)[0] for i in idx]

def initialize_pop(n: int, selfies_molecule, metric_function_list, mutation_function_list, generation=0):
    """
    Creates a data frame of n derived (SMILES) molecules, with evaluated metrics as additional columns
    """
    column_names = ["SMILES molecule", "Generation"] + [f"Metric {i+1}" for i in range(len(metric_function_list))]
    derivatives = pd.DataFrame(columns=column_names)

    while derivatives.drop_duplicates().shape[0] < n:
        molecules = [canonicalize_smiles(molecule) for molecule in
                     generate_derivatives(n - derivatives.drop_duplicates().shape[0], selfies_molecule, mutation_function_list)]
        
        data = [[molecule, generation] + [metric(molecule) for metric in metric_function_list] for molecule in molecules]

        derivatives = pd.concat([derivatives.drop_duplicates(), pd.DataFrame(data=data, columns=column_names)])
    
    return sort_df_by_col_index(derivatives, 1)

def populate_from_df(df: pd.DataFrame, n: int, metric_function_list: list, mutation_function_list: list, generation: int, 
                     fitness='Isolation score', remove_columns=['Isolation score','Pareto Front'], include_df=False): #TODO: implement `include df`
    """
    Outputs a populated data frame with molecules extracted from input data frame.
    """
    norm_weights = _normalize_weights(df, fitness).astype(int) 
    df.drop(remove_columns, axis=1, errors='ignore', inplace=True)
    seeds = get_smiles_list(df,verify=False) #seeds = staring molecules

    temp_gen_list = [initialize_pop(weighted_n, sf.encoder(seed), metric_function_list, mutation_function_list, generation) 
                     for seed, weighted_n in zip(seeds, norm_weights)]

    new_gen = pd.concat(temp_gen_list)
    return sort_df_by_col_index(new_gen, 1)

def min_max_normalize(df, column_name):
    df[column_name] = (df[column_name] - df[column_name].min()) / (df[column_name].max() - df[column_name].min())
    return df





def generate_derivatives_batch(n:int,
                               selfies_seeds: List[str],
                               mutation_function_list: list,
                               crossover = False,
                               unique = False,
                               crossover_ratio = 0.25,
                               crossover_type = 0) -> list:  #TODO implement `unique`
    """
    Creates `n` derivatives of seed SELFIES molecules, outputs SMILES list.
    If `crossover = True`, a `crossover_ratio` of `n` is taken and crossed over instead of mutated
    """

    

    if crossover:
        random.shuffle(selfies_seeds)
        n_crossed = round(n * crossover_ratio)
        seeds_crossed = selfies_seeds[:n_crossed]
        seeds_crossed = (xo.random_individual_crossover(seeds_crossed,n_crossed, crossover_type=crossover_type)) 
        n_untouched = n - n_crossed
        seeds_untouched = selfies_seeds[n_crossed:]
        output = selfies_to_smiles_list([shorten_rings(s) for s in seeds_crossed])
    else:
        seeds_untouched = selfies_seeds
        output = []



    #print("after:\n",seeds)
    len_seeds_untouched = len(seeds_untouched) #10
    partitions = len(mutation_function_list) #3

    idx = np.arange(len_seeds_untouched) % partitions #[0 1 2 0 1 2 0 1 2 0]


    for seed_i in range(len_seeds_untouched):
        i = idx[seed_i] # 0 or 1 or 2
        #TODO: implement `uniqueness` (non-repeating) (WHILE NEW MOLECULE IN OUTPUT GENERATE ANOTHER ONE)
        #TODO delete these try excepts... #maybe don't lol, good for debugging
        try:
            a =(seeds_untouched[seed_i])
        except:
            print(f"len_seeds_untouched[seed_i] err\nseeds_untouched: {type(seeds_untouched)}\nseed_i: {type(seed_i)}")
        try:
            a = (mutation_function_list[i])
        except:
            print("mutation_function_list[i] err")
        try:
            a=(mutation_function_list[i](seeds_untouched[seed_i])[0])
        except:
            print("mutation_function_list[i](seeds[seed_i])[0] err")
        new = shorten_rings(mutation_function_list[i](seeds_untouched[seed_i])[1])
        output.append(sf.decoder(new)) #use function 0/1/2 on every 3rd seed to get smiles


    return output

def initialize_pop_batch(n: int, selfies_molecule, metric_function_list, mutation_function_list, generation=0):
    """
    Creates a data frame of n derived (SMILES) molecules, with evaluated metrics as additional columns
    """
    seeds = [selfies_molecule]
    derivatives = generate_derivatives_batch(n,seeds,mutation_function_list, crossover=False) #crossover MUST be false, this is only for generation from a single lead molecule
    return _common_loop_1(derivatives, metric_function_list,generation)


def populate_from_df_batch(df: pd.DataFrame,
                           n: int,
                           metric_function_list: list,
                           mutation_function_list: list,
                           generation: int,
                           fitness='Isolation score',
                           remove_columns=['Isolation score','Pareto Front'],
                           include_seeds=True,
                           crossover=True,
                           unique = False,
                           crossover_ratio = 0.25,
                           crossover_type = 0, fitness_proportional = True, get_cost = True):
    """
    Outputs a populated data frame with molecules extracted from input data frame.
    """
    seeds = smiles_to_selfies_list([randomize_smiles(s) for s in get_smiles_list(df,verify=False)]) #seeds = staring molecules

    # I WANT n OF THEM
    # HOW MANY MORE DO I NEED TO MAKE?
    # n MINUS WHAT I HAVE, len(seeds)
    # THAT WILL BE m

    if include_seeds:
        m = n - len(seeds)
    else:
        m = n

    if fitness_proportional:
        df = min_max_normalize(df, fitness) 
        df[fitness]+=0.1 #give non-0 chance to 'worst' solutions
        proportional_seeds = random.choices(seeds, weights=(df[fitness]**2).tolist(), k=m) #create fitness proportional list
    else:
        proportional_seeds = random.choices(seeds, k=m)
    derivatives = generate_derivatives_batch(m,proportional_seeds,
                                             mutation_function_list,
                                             crossover=crossover,
                                             unique = unique, 
                                             crossover_ratio = crossover_ratio,
                                             crossover_type = crossover_type)

    if include_seeds:
        derivatives.extend(selfies_to_smiles_list(seeds))


    #derivatives = list (of SMILES)
    return _common_loop_1(derivatives, metric_function_list,generation)


def convert_seeds_to_df(df: pd.DataFrame, metric_function_list: list):
    """
    Takes df of seeds and adds metrics.
    """
    return _common_loop_1(get_smiles_list(df),metric_function_list,0)[0]

evaluated_canonical_chems = {} #TODO store chemicals with their metric values e.g. {'COC': [0.754, 0.598]}

from collections import Counter
def _common_loop_1(SMILES: List[str],metric_function_list: List,generation: int): #TODO: rename
    """
    Outputs SMILES 
    """


    canonical_smiles = [canonicalize_smiles(s) for s in SMILES] # 'C1C(C)C1', 'C1(C)CC1' -> if same, turn into unified form 
    chemical_counts = count_list_elements(canonical_smiles)     # {'C1C(C)C1': 5, 'COC': 3 , ...}
    chemical_counts_keys = list(chemical_counts.keys())         # ['C1C(C)C1', 'COC']
    chemical_counts_values = list(chemical_counts.values())     # [5, 3]
    metrics = []
    for metric in metric_function_list:
        try:
            metrics.append(metric(chemical_counts_keys)) #['C1CC1', 'COC'] -> [0.546, 0.895]
        except:
            print(f"metric(SMILES) err for:\nmetric: {metric}\nSMILES {chemical_counts_keys}")
    smiles_list, metric_list = scramble_quantified_smiles(chemical_counts, metrics[0]) #unified forms -> [ [random froms * n], [metrics * n] ]
    generations = [generation] * len(smiles_list) #list (of numeric values)
    data = [smiles_list, metric_list]#.copy
    data.insert(1,generations) # [[random froms * n], [metrics * n]] -> [[random froms * n], [generations], [metrics * n]]
    column_names = ["SMILES molecule", "Generation"] + [f"Metric {i+1}" for i in range(len(metric_function_list))]
    cost = len(metrics[0])
    
    return _lists_to_dataframe(data, column_names), cost

def count_list_elements(lst):
    return dict(Counter(lst))

def scramble_quantified_smiles(smiles_frequencies, metrics):
    smiles_list = []
    metric_list = []
    for (k, v), metric_val in zip(smiles_frequencies.items(), metrics):
        smiles_list.extend([k]*v)
        metric_list.extend([metric_val]*v)
    smiles_list = [randomize_smiles(s) for s in smiles_list]
    
    return smiles_list, metric_list




def _add_metrics(SMILES: List[str],metric_function_list: List): #TODO: remake wrt above function
    """
    Outputs SMILES 
    """
    metrics = []
    for metric in metric_function_list:
        try:
            metrics.append(metric(SMILES))
        except:
            print(f"metric(SMILES) err for:\nmetric: {metric}\nSMILES: {SMILES}")
    
    data = [SMILES]
    data.extend(metrics)
    column_names = ["SMILES molecule"] + [f"Metric {i+1}" for i in range(len(metric_function_list))]
    
    return _lists_to_dataframe(data, column_names)




def _lists_to_dataframe(data, column_names):
    # Check if data and column_names have same length
    try:
        if len(data) != len(column_names):
            raise ValueError('Number of columns in data does not match number of column names')
    except:
        print(f"{type(data)}, {type(column_names)}")
    data_dict = {column_name: column_data for column_name, column_data in zip(column_names, data)}
    return pd.DataFrame(data_dict)


def smiles_to_selfies_list(smiles: list) -> list:
    output = []
    for s in smiles:
        try:
            output.append(sf.encoder(s))
        except:
            print(f"Error converting SMILES {s} to SELFIES, so it was skipped")
            continue
    return output
    #return [sf.encoder(s) for s in smiles]

def safe_encoder(smiles):
    try:
        return sf.encoder(smiles)
    except Exception:
        return None

def process_df(df, column_name):
    with Pool(processes=4) as pool:  # adjust the number of processes based on your CPU
        df['encoded'] = pool.map(safe_encoder, df[column_name])

    df = df.dropna(subset=['encoded'])
    df = df.drop(columns=['encoded'])
    
    return df



def selfies_to_smiles_list(smiles: list) -> list: #veeery slow
    return [sf.decoder(s) for s in smiles]


def _normalize_weights(df: pd.DataFrame, weights='Isolation score', round = True):
    if round:
        return np.round(np.divide(df[weights]*100, df[weights].sum()), 0).astype(int)
    return np.divide(df[weights]*100, df[weights].sum())
def _dominates(row1, row2, minimize=True):
    if minimize:
        return np.all(row1 <= row2) and np.any(row1 < row2)
    else:
        return np.all(row1 >= row2) and np.any(row1 > row2)
    

def _classify_row_fronts(df, metrics, front, minimize=True):#TODO: change 'minimize' to List[bool]
    not_classified = df[df['Pareto Front'] == PLACEHOLDER_VALUE]
    for idx in not_classified.index:
        row = not_classified.loc[idx, metrics].values
        other_rows = not_classified.loc[not_classified.index != idx, metrics].values
        dominated_by_any = np.any(np.apply_along_axis(_dominates, 1, other_rows, row, minimize))
        if not dominated_by_any:
            df.at[idx, 'Pareto Front'] = front
    return df


def get_pareto_optimal(df, column_names, minimize=True):#TODO: change 'minimize' to List[bool]
    """
    Extracts Pareto optimal solutions from data frame
    """
    df['Pareto Front'] = PLACEHOLDER_VALUE
    df = _classify_row_fronts(df, column_names, 1, minimize)
    return df.loc[df['Pareto Front'] == 1]


def classify_all_pareto_fronts(df, column_names, minimize=True):#TODO: change 'minimize' to List[bool]
    """
    Classifies all solutions in data frame into Pareto fronts
    """
    df['Pareto Front'] = PLACEHOLDER_VALUE
    front = 1
    while df[df['Pareto Front'] == PLACEHOLDER_VALUE].shape[0] > 0:
        df = _classify_row_fronts(df, column_names, front, minimize)
        front += 1
    return df


def _get_euclidean_distance(row1, row2):
    return np.linalg.norm(row1 - row2)



def get_percent_best(df: pd.DataFrame, column_name, percentage, minimize=True): #FIXME works only for single column #TODO make a multi-column version
    if len(column_name)>1:
        print("too many columns to optimize, pick only one")
    return df.sort_values(column_name,ascending=minimize).head(int(np.ceil(len(df)*(percentage))))
    
def get_canonical_percent_best(df: pd.DataFrame, column_name, percentage, minimize=True): #FIXME works only for single column #TODO make a multi-column version
    if len(column_name)>1:
        print("too many columns to optimize, pick only one")
    df['SMILES molecule'] = (df['SMILES molecule'].apply(lambda x: canonicalize_smiles(x)))
    df = df.drop_duplicates()
    return df.sort_values(column_name,ascending=minimize).head(int(np.ceil(len(df)*(percentage))))

def get_percent_worst(df: pd.DataFrame, column_name, percentage, minimize=True): #FIXME works only for single column #TODO make a multi-column version
    if len(column_name)>1:
        print("too many columns to optimize, pick only one")
    return df.sort_values(column_name,ascending=minimize).tail(int(np.ceil(len(df)*(percentage))))

def get_isolation(df: pd.DataFrame, metrics: List[str], alpha = 1):
    """
    Gets the isolation score of Pareto front classified solutions 
    (an inverse of crowdidng distance for easier calculations)
    """
    col_name = 'Isolation score'
    df.insert(1, col_name, PLACEHOLDER_VALUE)

    while df[df[col_name] == PLACEHOLDER_VALUE].shape[0] > 0:
        not_classified = df[df[col_name] == PLACEHOLDER_VALUE].sort_index()
        index_list = not_classified.index.tolist()
        for idx in range(len(index_list)):
            if (index_list[idx] == df.index.min()) or (index_list[idx] == df.index.max()):
                df.at[index_list[idx], col_name] = 0.00001
            else:
                dist1 = _get_euclidean_distance(np.asarray(df.loc[index_list[idx], metrics].values), np.asarray(df.loc[index_list[idx-1], metrics].values))
                dist2 = _get_euclidean_distance(np.asarray(df.loc[index_list[idx], metrics].values), np.asarray(df.loc[index_list[idx+1], metrics].values))
                df.at[index_list[idx], col_name] =  ((dist1+dist2))
    df.at[df.index.min(), col_name] = df.at[df.index.max(), col_name] = df[col_name].max() * alpha # to ensure higher/lower chance of being picked
    df = df.sort_values(metrics[0])
    return df


def distribute_budget(b, i):
    # Create an array of weights in reverse order
    weights = list(range(i, 0, -1))

    # Normalize the weights so they sum up to 1
    norm_weights = [float(w) / sum(weights) for w in weights]

    # Calculate number of evaluations for each iteration based on the weights and the total budget
    num_evals = [int(b * w) for w in norm_weights]

    # Adjust the allocations to ensure the total budget is not exceeded due to rounding
    while sum(num_evals) > b:
        max_index = num_evals.index(max(num_evals))
        num_evals[max_index] -= 1

    # Check if we underused the budget due to rounding down and if so, add the difference to the first iteration
    if sum(num_evals) < b:
        num_evals[0] += b - sum(num_evals)

    return num_evals

def get_diversity(df: pd.DataFrame, generation:int):
    data = df[df['Generation']==generation]['SMILES molecule']
    return len(pd.unique(data.apply(lambda x: canonicalize_smiles(x))))/len(data)
    

def get_last_diversity(df: pd.DataFrame):
    data = df[df['Generation']==df['Generation'].max()]['SMILES molecule']
    return (len(pd.unique(data.apply(lambda x: canonicalize_smiles(x))))/len((data)))

import rdkit


#ADAPTED FROM STONED PAPER
def randomize_smiles(smiles_string):
    # convert SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles_string)
    '''Returns a random (dearomatized) SMILES given an rdkit mol object of a molecule.
    Parameters:
    mol (rdkit.Chem.rdchem.Mol) :  RdKit mol object (None if invalid smile string smi)
    
    Returns:
    smiles (string) : a SMILES string
    '''
    Chem.Kekulize(mol)
    return rdkit.Chem.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False,  kekuleSmiles=True) 

def shorten_rings(selfies):
    return selfies
    split_selfies_list = list(sf.split_selfies(selfies))
    ring_list = ['[Ring1]', '[Ring2]', '[Branch1]', '[=Branch1]']
    for i in range(len(split_selfies_list) - 1): # -1 to avoid index out of range error
        if split_selfies_list[i] == '[Ring1]' and split_selfies_list[i + 1] not in ring_list:
            split_selfies_list[i + 1] = '[=Branch1]'
    return ''.join(split_selfies_list)


