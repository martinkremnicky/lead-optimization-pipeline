import selfies as sf
import numpy as np
import pandas as pd
from rdkit import Chem
from typing import Tuple

ALPHABET = sf.get_semantic_robust_alphabet() # Gets the alphabet of robust symbols
PLACEHOLDER_VALUE = -2

def canonicalize_smiles(smiles):
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)
    except:
        print(f"Failed at {smiles}")
def canonicalize_selfies(selfies: str):
    return Chem.MolToSmiles(Chem.MolFromSmiles(sf.decoder(selfies)), canonical=True)

def validate(selfies_molecule:str)->Tuple[str, str]:
    """
    Converts a (generated) SELFIES molecule to SMILES and back to 'validate' it
    ### Output:
    tuple of (SMILES,SELFIES) notations of the same molecule
    """
    conversion_smi = sf.decoder(selfies_molecule)
    conversion_sf = sf.encoder(conversion_smi)
    return conversion_smi, conversion_sf

def generate_derivatives(n:int, selfies_molecule: str, mutation_function_list: list) -> list:
    """
    Creates derivatives of SELFIES molecule, outputs SMILES list by default.
    """
    idx = np.arange(n) % len(mutation_function_list)
    return [mutation_function_list[i](selfies_molecule)[0] for i in idx]

def sort_df_by_col_index(df: pd.DataFrame, index: int) -> pd.DataFrame:
    return df.sort_values(df.columns[index]).reset_index(drop=True)

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
                     weights='Isolation score', remove_columns=['Isolation score','Pareto Front'], include_df=False):
    """
    Outputs a populated data frame with molecules extracted from input data frame.
    """
    norm_weights = _normalize_weights(df, weights).astype(int)
    df.drop(remove_columns, axis=1, errors='ignore', inplace=True)
    molecules = df['SMILES molecule']

    temp_gen_list = [initialize_pop(weighted_n, sf.encoder(molecule), metric_function_list, mutation_function_list, generation) 
                     for molecule, weighted_n in zip(molecules, norm_weights)]

    new_gen = pd.concat(temp_gen_list)
    return sort_df_by_col_index(new_gen, 1)

def _normalize_weights(df: pd.DataFrame, weights='Isolation score'):
    return np.round(np.divide(df[weights]*100, df[weights].sum()), 0).astype(int)

def _dominates(row1, row2, minimize=True):
    if minimize:
        return np.all(row1 <= row2) and np.any(row1 < row2)
    else:
        return np.all(row1 >= row2) and np.any(row1 > row2)

def _classify_row_fronts(df, metrics, front, minimize=True):
    not_classified = df[df['Pareto Front'] == PLACEHOLDER_VALUE]
    for idx in not_classified.index:
        row = not_classified.loc[idx, metrics].values
        other_rows = not_classified.loc[not_classified.index != idx, metrics].values
        dominated_by_any = np.any(np.apply_along_axis(_dominates, 1, other_rows, row, minimize))
        if not dominated_by_any:
            df.at[idx, 'Pareto Front'] = front
    return df

def get_pareto_optimal(df, metrics, minimize=True):
    """
    Extracts Pareto optimal solutions from data frame
    """
    df['Pareto Front'] = PLACEHOLDER_VALUE
    df = _classify_row_fronts(df, metrics, 1, minimize)
    return df.loc[df['Pareto Front'] == 1]

def classify_all_pareto_fronts(df, metrics, minimize=True):
    """
    Classifies all solutions in data frame into Pareto fronts
    """
    df['Pareto Front'] = PLACEHOLDER_VALUE
    front = 1
    while df[df['Pareto Front'] == PLACEHOLDER_VALUE].shape[0] > 0:
        df = _classify_row_fronts(df, metrics, front, minimize)
        front += 1
    return df

def _get_euclidean_distance(row1, row2):
    return np.linalg.norm(row1 - row2)

def get_isolation(df: pd.DataFrame, metrics):
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
                dist1 = _get_euclidean_distance(df.loc[index_list[idx], metrics].values, df.loc[index_list[idx-1], metrics].values)
                dist2 = _get_euclidean_distance(df.loc[index_list[idx], metrics].values, df.loc[index_list[idx+1], metrics].values)
                df.at[index_list[idx], col_name] =  ((dist1+dist2))
    df.at[df.index.min(), col_name] = df.at[df.index.max(), col_name] = df[col_name].max() # * alpha # to ensure higher/lower chance of being picked
    df = df.sort_values(metrics[0])
    return df
