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
import copy
import numpy as np
from typing import List, Tuple, Dict
from collections import Counter
from copy import deepcopy
from pandas import DataFrame #FIXME not needed
import itertools
from collections import Counter
import edlib
#from numba import njit



def canonicalize_smile(smile: str) -> str:
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smile), canonical=True)
    except:
        print(f"Failed to 'canonicalize_smile' for {smile}")
        pass

def canonicalize_selfie(selfie: str) -> str:
    """returns SMILES"""
    try: 
        return Chem.MolToSmiles(Chem.MolFromSmiles(sf.decoder(selfie)), canonical=True)
    except:
        print(f"Failed to 'canonicalize_selfie' for {selfie}")
        pass

def validate(selfie:str)->Tuple[str, str]:
    """
    Converts a (generated) SELFIES molecule to SMILES and back to 'validate' it
    ### Output:
    tuple of (SMILES,SELFIES) notations of the same molecule
    """
    conversion_smi = sf.decoder(selfie)
    try:
        conversion_sf = sf.encoder(conversion_smi)
        return conversion_smi, conversion_sf
    except:
        print(f"functions.py:validate(selfie) - conversion_sf of smiles ({conversion_smi}) resulted in err (ORIGINAL SELFIES: {selfie})")

def sort_df_by_col_index(df: pd.DataFrame, index: int) -> pd.DataFrame:
    return df.sort_values(df.columns[index]).reset_index(drop=True)

def column_to_list(df: pd.DataFrame, col_name:str) -> pd.DataFrame:
    if col_name not in df.columns:
        raise ValueError(f"The DataFrame does not contain a column named '{col_name}'.")
    return df[col_name].tolist()

def get_smiles(df: pd.DataFrame, verify=True) -> list:
    accepted_column_names = ["SMILES","smiles", "SMILES"] #FIXME perhaps a more flexible approach would be better
    for name in accepted_column_names:
        try:
            if verify:
                return column_to_list(process_df(df,name),name)
            else:
                return column_to_list(df,name)
        except:
            pass
    print(f"Could not find a column with an accepted name. Accepted names:\n{accepted_column_names}") #TODO make `raise`?

def min_max_normalize(series: pd.Series, bias = 0) -> pd.Series:
    """
    Perform min-max normalization on a specific column in a DataFrame.

    Args:
        series (pd.Series): The input Series containing the data.
        minimal_value (float): The minimal value to add after normalization.

    Returns:
        pd.Series: The normalized values of the specified column, with the minimal value added.
    """
    series = np.asarray(series)
    #print(type(series),np.min(series))
    #print(f"{df[column_name][0]} {type(df[column_name][0])}\n{df[column_name].max()} {type(df[column_name].max())}\n{df[column_name].min()} {type(df[column_name].min())}")
    divisor = (np.max(series) - np.min(series))
    if divisor == 0:
        divisor = 1
    return np.divide((series - np.min(series)), divisor) + bias #give non-0 chance to 'worst' solutions



def generate_derivatives(n:int,
                        selfies_seeds: List[str],
                        mutations: list,
                        crossover = False,
                        crossover_ratio = 0.25,
                        crossover_type = 0,
                        short_rings=False) -> List[str]:  
    """
    Creates `n` derivatives of seed SELFIES molecules, outputs SMILES list.
    If `crossover = True`, a `crossover_ratio` of `n` is taken and crossed over instead of mutated
    """
    if True:
        seeds_untouched = []
        seeds_crossed = []
        output = []

        if crossover_ratio > 0:
            random.shuffle(selfies_seeds)
            n_cross = round(n * crossover_ratio)

            if n_cross > 1:
                seeds_to_cross = selfies_seeds[:n_cross]
                seeds_crossed = xo.random_individual_crossover(seeds_to_cross, n_cross, crossover_type=crossover_type)
                seeds_untouched = selfies_seeds[n_cross:]
            else:
                seeds_untouched = selfies_seeds
            output = selfies_to_smiles_list([shorten_rings(s) if short_rings else s for s in seeds_crossed])
        else:
            seeds_untouched = selfies_seeds


    #cycle = itertools.cycle(range(len(mutations)))
    len_su = len(seeds_untouched)
    mut_indexes = random.choices(range(len(mutations)),k=len_su)
    for i in range(len_su):
    #for seed in seeds_untouched:
        #i = next(cycle)
        #new = mutations[i](seed)[1]
        new = mutations[mut_indexes[i]](seeds_untouched[i])[1]
        output.append(sf.decoder(new))
    return output 

def initialize_pop(n: int, selfie, metrics, mutations, generation=0):
    """
    Creates a data frame of n derived (SMILES) molecules, with evaluated metrics as additional columns
    """
    seeds = [selfie]
    derivatives = generate_derivatives(n,seeds,mutations, crossover=False) #crossover MUST be false, this is only for generation from a single lead molecule #..does it?
    metric_evaluated_df, cost, lookup = evaluate_smiles(derivatives, metrics,generation)
    fitness_evaluated_df = add_fitness(metric_evaluated_df, generation)
    return fitness_evaluated_df, cost, lookup
    
def populate_from_df(df: pd.DataFrame,
                    n: int,
                    metrics: list,
                    mutations: list,
                    generation: int,
                    fitness='Fitness',
                    include_seeds=True,
                    crossover=True,
                    crossover_ratio = 0.25,
                    crossover_type = 0,
                    proportional_sampling = True,
                    lookup_dict = {},
                    randomize_seeds = False,#True,
                    delete_duplicates = False,
                    exp = 2,
                    min_fitness = 0.1,
                    normalize_fitness = True,
                    avg_dist = True,
                    use_fitness = True):
    """
    Outputs a populated data frame with molecules extracted from input data frame.
    """
    #print('fn.py - populate_from_df',metrics)
    df.sort_values(fitness,ascending=False, inplace=True)

    if randomize_seeds:
        df["SMILES"] = df["SMILES"].apply(randomize_smile)
    #    seeds = smiles_to_selfies_list(get_smiles(df,verify=False)) #seeds = starting molecules
    else:
        if delete_duplicates:
            df.drop_duplicates('Canonical SMILES',inplace=True)
    seeds = smiles_to_selfies_list(get_smiles(df,verify=False))
    
    
    # I WANT n OF THEM
    # HOW MANY MORE DO I NEED TO MAKE?
    # n MINUS WHAT I HAVE, len(seeds)
    # THAT WILL BE m

    min_n = 2 #has to be above 1, else infinite loop due to elitism

    if n<min_n:
        n = min_n

    if include_seeds:
        m = n - len(seeds)
    else:
        m = n
    #print('seeds:',len(seeds),m,'n',n,'generation',generation,)

    if normalize_fitness:
        df[fitness] = min_max_normalize(df[fitness], min_fitness)

    if exp != 1:
        df[fitness]=np.float_power(df[fitness],np.float16(exp))


    if proportional_sampling:
        sampled_seeds = random.choices(seeds, weights=(df[fitness]), k=m) #create fitness proportional list
    else:
        sampled_seeds = random.choices(seeds, k=m)
    derivatives = generate_derivatives(n=m,
                                        selfies_seeds= sampled_seeds,
                                        mutations=mutations,
                                        crossover=crossover,
                                        crossover_ratio = crossover_ratio,
                                        crossover_type = crossover_type)
    if include_seeds:
        derivatives.extend(selfies_to_smiles_list(seeds))

    metric_evaluated_df, cost, lookup = evaluate_smiles(derivatives, metrics, generation, lookup_dict, randomize_seeds)
    if use_fitness:
        fitness_evaluated_df = add_fitness(metric_evaluated_df, avg_dist, generation)
        #print(f"+=+=+ der, met, fit lens: {len(derivatives),len(metric_evaluated_df),len(fitness_evaluated_df)}")
        return fitness_evaluated_df, cost, lookup
    return metric_evaluated_df, cost, lookup

def seeds_to_pop(df: pd.DataFrame, metrics: List):
    """
    Takes df of seeds and adds metrics as columns.
    """
    #print('fn.py - seeds_to_pop',metrics)
    metric_evaluated_df, cost, lookup = evaluate_smiles(get_smiles(df),metrics,0)
    fitness_evaluated_df = add_fitness(metric_evaluated_df, True, 0)
    return fitness_evaluated_df, cost, lookup

def dataframe_to_dict(df): #TODO: un-hardcode
    return df.set_index('Canonical SMILES')['Metric 1'].to_dict()


def evaluate_smiles(smiles: List[str], metrics: List, generation: int, lookup = {}, randomize =True): #TODO: rename #FIXME what if more than one metric is presented?
    #print('fn.py - evaluate_smiles',metrics)
    """
    Outputs:
    1. SMILES df with metrics
    2. evaluation costs
    3. updated metric lookup
    """
    canonical_smiles = [canonicalize_smile(s) for s in smiles]  # 'C1C(C)C1', 'C1(C)CC1' -> if same, turn into unified form 
    chemical_freq = Counter(canonical_smiles)     # {'C1C(C)C1': 5, 'COC': 3 , ...}

    new_chems = []
    old_chems = []

    [old_chems.append(c) if c in lookup.keys() else new_chems.append(c) for c in chemical_freq.keys() ]
    
    old_metric_values = [lookup[c] for c in old_chems]
    new_metric_values = []
    for metric in metrics:
        try:
            new_metric_values.append(metric(new_chems))  #['C1CC1', 'COC'] -> [0.546, 0.895]
        except:
            print(f"metric(SMILES) err for:\nmetric: {metric}\nSMILES {chemical_freq.keys()}")
            new_metric_values.append(metric(new_chems))
    
    old_chem_freq = {k: v for k, v in chemical_freq.items() if k in old_chems}
    new_chem_freq = {k: v for k, v in chemical_freq.items() if k in new_chems}

    merged_freq = {**new_chem_freq, **old_chem_freq}
    merged_metrics = new_metric_values[0] + old_metric_values #FIXME why the [0]? is it bcs it's a single metric? then it needs fixing

    smiles_list, metric_list = expand_by_frequency(merged_freq, merged_metrics, randomize) #unified forms -> [ [random froms * n], [metrics * n] ]
    selfies_list = [sf.encoder(sm) for sm in smiles_list]
    can_smiles = [canonicalize_smile(sm) for sm in smiles_list]
    can_selfies = [sf.encoder(canonicalize_smile(sm)) for sm in smiles_list]
    
    split_selfies = [list(sf.split_selfies(sel)) for sel in can_selfies]

    generations = [generation] * len(smiles_list) 
    data = [smiles_list, selfies_list, can_smiles, can_selfies, split_selfies, generations, metric_list]
    
    column_names = ["SMILES","SELFIES","Canonical SMILES","Canonical SELFIES", "Split canonical SELFIES", "Generation"] + [f"Metric {i+1}" for i in range(len(new_metric_values))]
    
    cost = len(new_metric_values[0]) #FIXME the [0]..

    lookup.update({nc: nmv for nc, nmv in zip(new_chems, new_metric_values[0])})


    metric_evaluated_df = _lists_to_dataframe(data, column_names)
    return metric_evaluated_df, cost, lookup


def add_fitness(df:pd.DataFrame, avg_dist = True, generation = 1): #FIXME hardcoded for Metric 1
    
    
    c = 25
    diff_factor = 0.1
    low_bar = 1.0

    if (not avg_dist) or (generation>c):
        df['Fitness'] = df['Metric 1']
        df['avg_distance'] = 0
        df['time_decline'] = 0
        df['time_increase'] = 0
        df['rep'] = 0
        df['atr'] = 0
        df['add'] = 0
    else:
        repulsion = 0
        attraction = 0
        time_decline = (-(np.divide(1,c)*generation)+1)
        #time_decline = 1
        time_increase = generation-c




        column_name = "Split canonical SELFIES"
        metric_ones = df['Metric 1'].values
        metric_max = metric_ones.max()
        max_addition = metric_max  * time_decline * diff_factor
        above_threshold_selfies = df[(df['Metric 1']+max_addition)>metric_max * low_bar][column_name]
        above_threshold_selfies = [tuple(s) for s in above_threshold_selfies]



        if False:
            SELFIES_as_arr = df[column_name].values
            avg_distance = df[column_name].apply(avg_edit_distance, df_as_arr=SELFIES_as_arr)
            avg_distance = min_max_normalize(avg_distance) #arr of 0 to 1
            avg_distance = np.multiply((avg_distance),np.float16(diff_factor))
        else:
            SELFIES_as_arr = df[column_name].values
            avg_distance = df[column_name].apply(
                lambda selfie: (avg_edit_distance(row= selfie,df_as_arr=SELFIES_as_arr)) if tuple(selfie) in above_threshold_selfies else 0


            )
            avg_distance = min_max_normalize(avg_distance) #arr of 0 to 1
            avg_distance = np.multiply((avg_distance),np.float16(diff_factor)) #FUUUUUUUUUUUCK
        repulsion = avg_distance**2 * time_decline
        
        
        
        # repulsion =  (np.multiply((  min_max_normalize(avg_distance)  ),np.float16(diff_factor)))**2 * time_decline
        
        
        
        #if generation>25:
        #    print('A:\n',avg_distance)
        #total_distance = df[column_name].apply(total_edit_distance, df_as_arr=df_as_arr)
        #total_distance = min_max_normalize(total_distance.to_numpy())
        #total_distance = (total_distance.to_numpy())*0.01
        #avg_distance = np.divide(total_distance, np.size(total_distance))
        #avg_distance = np.mean(total_distance)

        #if generation>25:
        #    print('B:\n',avg_distance)

        #if generation>25:
        #    print('C:\n',avg_distance)
        #repulsion = ((100/(generation+1)) * avg_distance) 
        #attraction = (1.05**generation * avg_distance) *0.1

        k1 = 1
        k2 = 2
        k3 = 3
        k4 = 4
        




        add = df['Metric 1'].max() * np.asarray([r if r>0 else 0 for r in repulsion])
        df['avg_distance'] = avg_distance
        df['time_decline'] = time_decline
        df['time_increase'] = time_increase
        df['rep'] = repulsion
        df['atr'] = attraction
        df['add'] = add




        df['Fitness'] = df['Metric 1'] + add
    return df.sort_values('Fitness',ascending=False)

def avg_edit_distance(row, df_as_arr):
    similarity_sum = 0
    num_comparisons = 0
    for other_string in df_as_arr:
        if row != other_string:
            similarity_sum += edlib.align(((row)),((other_string)))['editDistance']
            num_comparisons += 1
    return similarity_sum / num_comparisons if num_comparisons > 0 else 0

#def avg_edit_distance(row, df_as_arr):
#    distances = [edlib.align(row, other_string)['editDistance'] for other_string in df_as_arr if row != other_string]
#    return np.mean(distances) if distances else 0



def total_edit_distance(row, df_as_arr):
    similarity_sum = sum(edlib.align(row, other_string)['editDistance'] for other_string in df_as_arr if row != other_string)
    return similarity_sum


import numpy as np
from itertools import combinations

def unique_handshakes(arr):
    # Convert the NumPy array to a Python set to remove duplicates
    unique_values = set(arr)
    
    # Create the unique handshakes using combinations
    handshakes = list(combinations(unique_values, 2))
    
    return handshakes




def count_list_elements(lst):
    return dict(Counter(lst))

def expand_by_frequency(smiles_frequencies, metrics: List, randomize = True):
    smiles_list = []
    metric_list = []
    for (k, v), metric_val in zip(smiles_frequencies.items(), metrics):
        smiles_list.extend([k]*v)
        metric_list.extend([metric_val]*v)
    if randomize:
        smiles_list = [randomize_smile(s) for s in smiles_list]
    
    return smiles_list, metric_list

def filter_dictionary(frequencies, selected_keys):
    selected_dict = {}
    for key in selected_keys:
        if key in frequencies:
            selected_dict[key] = frequencies[key]
    return selected_dict


def _add_metrics(smiles: List[str], metrics: List): #TODO: remake wrt above function
    """
    Outputs SMILES 
    """
    metrics = []
    for metric in metrics:
        try:
            metrics.append(metric(smiles))
        except:
            print(f"metric(SMILES) err for:\nmetric: {metric}\nSMILES: {smiles}")
    
    data = [smiles]
    data.extend(metrics)
    column_names = ["SMILES"] + [f"Metric {i+1}" for i in range(len(metrics))]
    
    return _lists_to_dataframe(data, column_names)

def _lists_to_dataframe(data: List[List], column_names: List[str]):
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
        except: #FIXME delet this
            print(f"Error converting SMILES {s} to SELFIES, so it was skipped")
            continue
    return output
    #return [sf.encoder(s) for s in smiles]

def safe_encoder(smiles):
    try:
        return sf.encoder(smiles)
    except Exception:
        return None

def process_df(df:pd.DataFrame, column_name: str) -> pd.DataFrame:
    with Pool(processes=4) as pool:  # adjust the number of processes based on your CPU
        df['encoded'] = pool.map(safe_encoder, df[column_name])

    df = df.dropna(subset=['encoded'])
    df = df.drop(columns=['encoded'])
    return df

def selfies_to_smiles_list(smiles: List[str]) -> List[str]: #veeery slow
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

def get_pareto_optimal(df: pd.DataFrame, column_names: List[str], minimize=True) -> pd.DataFrame:#TODO: change 'minimize' to List[bool]
    """
    Extracts Pareto optimal solutions from data frame
    """
    df['Pareto Front'] = PLACEHOLDER_VALUE
    df = _classify_row_fronts(df, column_names, 1, minimize)
    return df.loc[df['Pareto Front'] == 1]


def classify_all_pareto_fronts(df: pd.DataFrame, column_names: List[str], minimize=True) -> pd.DataFrame:#TODO: change 'minimize' to List[bool]
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

def get_percent_best(df: pd.DataFrame, column_name: str, percentage: float, minimize=True): #FIXME works only for single column #TODO make a multi-column version
    if len(column_name)>1:
        print("too many columns to optimize, pick only one")
    return df.sort_values(column_name,ascending=minimize).head(int(np.ceil(len(df)*(percentage))))
    
def get_canonical_percent_best(df: pd.DataFrame, column_name: str, percentage: float, minimize=True): #FIXME works only for single column #TODO make a multi-column version
    if len(column_name)>1:
        print("too many columns to optimize, pick only one")
    #df['SMILES'] = (df['SMILES'].apply(lambda x: canonicalize_smile(x)))
    df = df.drop_duplicates(subset='Canonical SMILES')
    return df.sort_values(column_name,ascending=minimize).head(int(np.ceil(len(df)*(percentage))))

def get_percent_worst(df: pd.DataFrame, column_name: str, percentage: float, minimize=True): #FIXME works only for single column #TODO make a multi-column version
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
                df.at[index_list[idx], col_name] = 0.00001 #hmm
            else:
                dist1 = _get_euclidean_distance(np.asarray(df.loc[index_list[idx], metrics].values), np.asarray(df.loc[index_list[idx-1], metrics].values))
                dist2 = _get_euclidean_distance(np.asarray(df.loc[index_list[idx], metrics].values), np.asarray(df.loc[index_list[idx+1], metrics].values))
                df.at[index_list[idx], col_name] =  ((dist1+dist2))
    df.at[df.index.min(), col_name] = df.at[df.index.max(), col_name] = df[col_name].max() * alpha # to ensure higher/lower chance of being picked
    df = df.sort_values(metrics[0])
    return df

def distribute_budget(budget: int, iterations:int) -> List[int]:
    # Create an array of weights in reverse order
    weights = list(range(iterations, 0, -1))

    # Normalize the weights so they sum up to 1
    norm_weights = [float(w) / sum(weights) for w in weights]

    # Calculate number of evaluations for each iteration based on the weights and the total budget
    num_evals = [int(budget * w) for w in norm_weights]

    # Adjust the allocations to ensure the total budget is not exceeded due to rounding
    while sum(num_evals) > budget:
        max_index = num_evals.index(max(num_evals))
        num_evals[max_index] -= 1

    # Check if we underused the budget due to rounding down and if so, add the difference to the first iteration
    if sum(num_evals) < budget:
        num_evals[0] += budget - sum(num_evals)

    return num_evals

def get_diversity(df: pd.DataFrame, generation:int) -> float:
    data = df[df['Generation']==generation]['SMILES']
    return len(pd.unique(data.apply(lambda x: canonicalize_smile(x))))/len(data)
    

def get_last_diversity(df: pd.DataFrame) -> float:
    data = df[df['Generation']==df['Generation'].max()]['SMILES']
    return (len(pd.unique(data.apply(lambda x: canonicalize_smile(x))))/len((data)))

import rdkit


#ADAPTED FROM STONED PAPER
def randomize_smile(smile:str) -> str:
    """returns SMILES"""
    # convert SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smile)
    '''Returns a random (dearomatized) SMILES given an rdkit mol object of a molecule.
    Parameters:
    mol (rdkit.Chem.rdchem.Mol) :  RdKit mol object (None if invalid smile string smi)
    
    Returns:
    smiles (string) : a SMILES string
    '''
    Chem.Kekulize(mol)
    return rdkit.Chem.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False,  kekuleSmiles=True) 

def randomize_selfie(selfie:str) -> str:
    return sf.encoder(randomize_smile(sf.decoder(selfie)))



def shorten_rings(selfie) -> str:
    split_selfies_list = list(sf.split_selfies(selfie))
    ring_list = ['[Ring1]', '[Ring2]', '[Branch1]', '[=Branch1]']
    for i in range(len(split_selfies_list) - 1): # -1 to avoid index out of range error
        if split_selfies_list[i] == '[Ring1]' and split_selfies_list[i + 1] not in ring_list:
            split_selfies_list[i + 1] = '[=Branch1]'
    return ''.join(split_selfies_list)


