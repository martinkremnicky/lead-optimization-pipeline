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
from pandas import DataFrame



def canonicalize_smile(smile: str) -> str:
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smile), canonical=True)
    except:
        print(f"Failed to 'canonicalize_smile' for {smile}")
        pass

def canonicalize_selfie(selfie: str) -> str:
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
    accepted_column_names = ["SMILES molecule","smiles", "SMILES"] #FIXME perhaps a more flexible approach would be better
    for name in accepted_column_names:
        try:
            if verify:
                return column_to_list(process_df(df,name),name)
            else:
                return column_to_list(df,name)
        except:
            pass
    print(f"Could not find a column with an accepted name. Accepted names:\n{accepted_column_names}") #TODO make `raise`?

""" def generate_derivatives(n:int, selfie: str, mutations: list) -> list:
    """
    #Creates `n` derivatives of SELFIES molecule, outputs SMILES list by default.
"""
    idx = np.arange(n) % len(mutations) #e.g. [0 1 2 3 4 5...] % 3 -> [0 1 2 0 1 2 ...] 
    return [mutations[i](selfie)[0] for i in idx]

def initialize_pop(n: int, selfie, metrics, mutations, generation=0):
    """
    #Creates a data frame of n derived (SMILES) molecules, with evaluated metrics as additional columns
"""
    column_names = ["SMILES molecule", "Generation"] + [f"Metric {i+1}" for i in range(len(metrics))]
    derivatives = pd.DataFrame(columns=column_names)

    while derivatives.drop_duplicates().shape[0] < n:
        molecules = [canonicalize_smile(molecule) for molecule in
                     generate_derivatives(n - derivatives.drop_duplicates().shape[0], selfie, mutations)]
        
        data = [[molecule, generation] + [metric(molecule) for metric in metrics] for molecule in molecules]

        derivatives = pd.concat([derivatives.drop_duplicates(), pd.DataFrame(data=data, columns=column_names)])
    
    return sort_df_by_col_index(derivatives, 1)

def populate_from_df(df: pd.DataFrame, n: int, metrics: list, mutations: list, generation: int, 
                     fitness='Isolation score', remove_columns=['Isolation score','Pareto Front'], include_df=False): #TODO: implement `include df`
    """
    #Outputs a populated data frame with molecules extracted from input data frame.
"""
    norm_weights = _normalize_weights(df, fitness).astype(int) 
    df.drop(remove_columns, axis=1, errors='ignore', inplace=True)
    seeds = get_smiles(df,verify=False) #seeds = staring molecules

    temp_gen_list = [initialize_pop(weighted_n, sf.encoder(seed), metrics, mutations, generation) 
                     for seed, weighted_n in zip(seeds, norm_weights)]

    new_gen = pd.concat(temp_gen_list)
    return sort_df_by_col_index(new_gen, 1)
 """
def min_max_normalize(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
    df[column_name] = (df[column_name] - df[column_name].min()) / (df[column_name].max() - df[column_name].min())
    return df

def generate_derivatives(n:int,
                        selfies_seeds: List[str],
                        mutations: list,
                        crossover = False,
                        force_unique = False, #TODO implement `force_unique`
                        crossover_ratio = 0.25,
                        crossover_type = 0,
                        short_rings=False) -> List[str]:  
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
        if short_rings:
            output = selfies_to_smiles_list([shorten_rings(s) for s in seeds_crossed])
        else:
            output = selfies_to_smiles_list(seeds_crossed)
    else:
        seeds_untouched = selfies_seeds
        output = []

    len_seeds_untouched = len(seeds_untouched) #10
    partitions = len(mutations) #3

    idx = np.arange(len_seeds_untouched) % partitions #[0 1 2 0 1 2 0 1 2 0] #TODO change to cycle

    for seed_i in range(len_seeds_untouched):
        i = idx[seed_i] # 0 or 1 or 2
        #TODO: implement "uniqueness" (non-repeating) (WHILE NEW MOLECULE IN OUTPUT GENERATE ANOTHER ONE)
        #TODO delete these try excepts... #maybe don't lol, good for debugging
        try:
            a =(seeds_untouched[seed_i])
        except:
            print(f"seeds_untouched[seed_i] err\nseeds_untouched: {type(seeds_untouched)}\n{seeds_untouched}\nseed_i: {type(seed_i)}\n{seed_i}")
        try:
            a = (mutations[i])
        except:
            print("mutations[i] err")
            print(f"mut: {mutations}\ni: {i}")
        try:
            a=(mutations[i](seeds_untouched[seed_i])[0])
        except:
            print("mutations[i](seeds_untouched[seed_i])[0] err")
            print(f"mut: {mutations}\ni: {i} \nseeds_untouched {type(seeds_untouched)}\n{seeds_untouched}\nseed_i: {type(seed_i)}\n{seed_i}")
        if short_rings:
            new = shorten_rings(mutations[i](seeds_untouched[seed_i])[1])
        else:
            new = mutations[i](seeds_untouched[seed_i])[1] #use function 0/1/2 on every 3rd seed to get smiles
        output.append(sf.decoder(new))  
    return output

def initialize_pop(n: int, selfie, metrics, mutations, generation=0):
    """
    Creates a data frame of n derived (SMILES) molecules, with evaluated metrics as additional columns
    """
    seeds = [selfie]
    derivatives = generate_derivatives(n,seeds,mutations, crossover=False) #crossover MUST be false, this is only for generation from a single lead molecule
    return _common_loop_1(derivatives, metrics,generation)

def populate_from_df(df: pd.DataFrame,
                    n: int,
                    metrics: list,
                    mutations: list,
                    generation: int,
                    fitness='Isolation score',
                    include_seeds=True,
                    crossover=True,
                    force_unique = False,
                    crossover_ratio = 0.25,
                    crossover_type = 0,
                    fitness_proportional = True,
                    lookup_dict = {},
                    get_cost = True):
    """
    Outputs a populated data frame with molecules extracted from input data frame.
    """
    seeds = smiles_to_selfies_list([randomize_smiles(s) for s in get_smiles(df,verify=False)]) #seeds = staring molecules

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
        df[fitness]+=0.075 #give non-0 chance to 'worst' solutions
        proportional_seeds = random.choices(seeds, weights=(df[fitness]**2).tolist(), k=m) #create fitness proportional list
    else:
        proportional_seeds = random.choices(seeds, k=m)
    derivatives = generate_derivatives(n=m,
                                        selfies_seeds= proportional_seeds,
                                        mutations=mutations,
                                        crossover=crossover,
                                        force_unique = force_unique, 
                                        crossover_ratio = crossover_ratio,
                                        crossover_type = crossover_type)
    if include_seeds:
        derivatives.extend(selfies_to_smiles_list(seeds))
    #derivatives = list (of SMILES)
    return _common_loop_1(derivatives, metrics, generation, lookup_dict)

def convert_seeds_to_df(df: pd.DataFrame, metrics: List):
    """
    Takes df of seeds and adds metrics.
    """
    return _common_loop_1(get_smiles(df),metrics,0)

evaluated_canonical_chems = {} #TODO store chemicals with their metric values e.g. {'COC': [0.754, 0.598]}

from collections import Counter
if True:
    def _common_loop_1(smiles: List[str], metrics: List, generation: int, lookup = {}): #TODO: rename #FIXME what if more than one metric is presented?
        """
        Outputs SMILES 
        """
        canonical_smiles = [canonicalize_smile(s) for s in smiles]  # 'C1C(C)C1', 'C1(C)CC1' -> if same, turn into unified form 
        chemical_freq = count_list_elements(canonical_smiles)     # {'C1C(C)C1': 5, 'COC': 3 , ...}
        chemical_freq_keys = list(chemical_freq.keys())         # ['C1C(C)C1', 'COC']
        #chemical_freq_values = list(chemical_counts.values())     # [5, 3]

        #print(f"lookup (len {len(lookup)}) {lookup}")
        new_chems = []
        old_chems = []
        old_metric_values = []
        new_metric_values = []
        for c in chemical_freq_keys:
            if c in lookup.keys():
                old_chems.append(c)
                old_metric_values.append(lookup[c])
            else:
                new_chems.append(c)

        for metric in metrics:
            try:
                new_metric_values.append(metric(new_chems))  #['C1CC1', 'COC'] -> [0.546, 0.895]
            except:
                print(f"metric(SMILES) err for:\nmetric: {metric}\nSMILES {chemical_freq_keys}")

        old_chem_freq = filter_dictionary(chemical_freq, old_chems)
        new_chem_freq = filter_dictionary(chemical_freq, new_chems)
        if old_chem_freq != {}:
            #print(len(old_chem_freq),'+',  len(new_chem_freq), '=')
            new_chem_freq.update(old_chem_freq)
            #print(len(new_chem_freq))
        merged_freq = copy.copy(new_chem_freq)
        #print(f"old_chem_freq\n{old_chem_freq}\nnew_chem_freq\n{new_chem_freq}\nmerged_freq\n{merged_freq}")
        merged_chems = new_chems + old_chems
        merged_metrics = new_metric_values[0] + old_metric_values #FIXME why the [0]? is it bcs it's a single metric? then it needs fixing
        

        smiles_list, metric_list = scramble_quantified_smiles(merged_freq, merged_metrics) #unified forms -> [ [random froms * n], [metrics * n] ]
        generations = [generation] * len(smiles_list) #list (of numeric values)
        data = [smiles_list, generations, metric_list]#.copy
        #data.insert(1,generations) # [[random froms * n], [metrics * n]] -> [[random froms * n], [generations], [metrics * n]]
        column_names = ["SMILES molecule", "Generation"] + [f"Metric {i+1}" for i in range(len(new_metric_values))]
        cost = len(new_metric_values[0]) #FIXME the [0]..

        #EXPAND LOOKUP DICT
        for nc, nmv in zip(new_chems,new_metric_values[0]): #FIXME the [0]
            lookup[nc] = nmv
        #print(f"new lookup (len {len(lookup)}) {lookup}")
        #print(f"cost {cost}")
        return _lists_to_dataframe(data, column_names), cost, lookup
else:

    def _common_loop_1(smiles: List[str], metrics: List, generation: int, lookup: Dict = None) -> Tuple[DataFrame, int, Dict]:
        """
        Outputs SMILES 
        """
        # If lookup is not provided, initialize it as an empty dictionary
        if lookup is None:
            lookup = {}

        canonical_smiles = [canonicalize_smile(s) for s in smiles]
        chemical_freq = Counter(canonical_smiles)

        old_chems = list(set(chemical_freq.keys()) & set(lookup.keys()))
        new_chems = list(set(chemical_freq.keys()) - set(old_chems))
        
        old_metric_values = [lookup[c] for c in old_chems]
        
        new_metric_values = []
        for metric in metrics:
            try:
                new_metric_values.append(metric(new_chems))
            except Exception as e:
                print(f"metric(SMILES) err for:\nmetric: {metric}\nSMILES {chemical_freq.keys()}")
                print(f"Error: {e}")
        
        old_chem_freq = {k: v for k, v in chemical_freq.items() if k in old_chems}
        new_chem_freq = {k: v for k, v in chemical_freq.items() if k in new_chems}

        merged_freq = {**new_chem_freq, **old_chem_freq}
        
        merged_chems = new_chems + old_chems
        merged_metrics = new_metric_values[0] + old_metric_values

        smiles_list, metric_list = scramble_quantified_smiles(merged_freq, merged_metrics)
        generations = [generation] * len(smiles_list)

        data = [smiles_list, generations, metric_list]
        column_names = ["SMILES molecule", "Generation"] + [f"Metric {i+1}" for i in range(len(new_metric_values))]

        cost = len(new_metric_values[0])
        
        #EXPAND LOOKUP DICT
        for nc, nmv in zip(new_chems, new_metric_values[0]):
            lookup[nc] = nmv

        return _lists_to_dataframe(data, column_names), cost, lookup


def count_list_elements(lst):
    return dict(Counter(lst))

def scramble_quantified_smiles(smiles_frequencies, metrics: List):
    smiles_list = []
    metric_list = []
    for (k, v), metric_val in zip(smiles_frequencies.items(), metrics):
        smiles_list.extend([k]*v)
        metric_list.extend([metric_val]*v)
    smiles_list = [randomize_smiles(s) for s in smiles_list]
    
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
    column_names = ["SMILES molecule"] + [f"Metric {i+1}" for i in range(len(metrics))]
    
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
    df['SMILES molecule'] = (df['SMILES molecule'].apply(lambda x: canonicalize_smile(x)))
    df = df.drop_duplicates()
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
    data = df[df['Generation']==generation]['SMILES molecule']
    return len(pd.unique(data.apply(lambda x: canonicalize_smile(x))))/len(data)
    

def get_last_diversity(df: pd.DataFrame) -> float:
    data = df[df['Generation']==df['Generation'].max()]['SMILES molecule']
    return (len(pd.unique(data.apply(lambda x: canonicalize_smile(x))))/len((data)))

import rdkit


#ADAPTED FROM STONED PAPER
def randomize_smiles(smile:str) -> str:
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
    return sf.encoder(randomize_smiles(sf.decoder(selfie)))



def shorten_rings(selfie) -> str:
    split_selfies_list = list(sf.split_selfies(selfie))
    ring_list = ['[Ring1]', '[Ring2]', '[Branch1]', '[=Branch1]']
    for i in range(len(split_selfies_list) - 1): # -1 to avoid index out of range error
        if split_selfies_list[i] == '[Ring1]' and split_selfies_list[i + 1] not in ring_list:
            split_selfies_list[i + 1] = '[=Branch1]'
    return ''.join(split_selfies_list)


