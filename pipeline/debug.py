from matplotlib import pyplot as plt
import selfies as sf
import mutations as mut
from functools import partial
import metrics as met
import functions as fn
from constants import *
from rdkit import Chem
from rdkit.Chem import Draw
import seaborn as sns
import pandas as pd
from datetime import datetime
import time
from tqdm import tqdm
from guacamol import standard_benchmarks
import numpy as np
import random
from copy import copy
import crossovers as xo

SEED_SIZE_LIST = [10, 100, 1000, 10000, 100000]
SEED_SIZE = SEED_SIZE_LIST[1]
SEED_LIST = [0,1,2]
SEED = SEED_LIST[2]
BUDGET_LIST = [100, 1000, 10000, 100000]
BUDGET = BUDGET_LIST[2]
GENERATIONS = 50




SAVE = True


#f = standard_benchmarks.zaleplon_with_other_formula().objective.score_list()
celecoxib = 'O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'
metric_function_list = [
    #partial(standard_benchmarks.hard_osimertinib().objective.score_list)
    partial(met.compound_similarity,target_smiles = celecoxib)
]
metrics = ['Metric 1']


mutation_function_list = [
    partial(mut.replacement)
    ,partial(mut.addition,fragment_size=1,rings=False)
    ,partial(mut.deletion,n=1)
]

current_dateTime = datetime.now()
time_format = f"{current_dateTime.date()}_{str(current_dateTime.hour).zfill(2)}-{str(current_dateTime.minute).zfill(2)}-{str(current_dateTime.second).zfill(2)}"

for SEED in [SEED_LIST[0]]:#SEED_LIST[1:]:
    for SEED_SIZE in [SEED_SIZE_LIST[2]]:#SEED_SIZE_LIST[1:]:#[SEED_SIZE_LIST[1]]:#
        SEED_PATH = f"./data/seed_{SEED}/rand_{SEED_SIZE}.tsv"
        seed_df = pd.read_table(SEED_PATH) 
        initial_pop = fn.convert_seeds_to_df(seed_df,metric_function_list)
        for BUDGET in BUDGET_LIST[2]:#BUDGET_LIST[:4]:#BUDGET_LIST[:4]#
            setup_name = f"seed_{SEED}_rand_{SEED_SIZE}"
            file_name = f"{time_format}_{setup_name}_budget_{BUDGET}"
            B0 = copy(BUDGET)
            print(f"Working with budget of {B0}")
            t0 = time.time()
            first_generation_fraction = 0.05 #hyperparam
            N = int(BUDGET * first_generation_fraction) #1st generation size

            if len(initial_pop)>=N:
                initial_best = initial_pop.head(N)
            else:
                initial_best = initial_pop
            gen_history = pd.DataFrame(initial_best)
            temp_best = initial_best.copy(deep=True)
            cost_history = []

            #rest of hyperparams
            crossover = True
            crossover_ratio = 0.2
            generation = 0
            
            next_generation_fraction = 0.05 #size of next generation; BUDGET * next_generation_fraction = new N
            initial_sample_fraction = 1.5 #initial_sample_size = len(temp_best) *  initial_sample_fraction ))
            minimal_next_generation_fraction = 0.0025
            max_gens = 500
            percent_best_fraction = 0.1
            include_initial_pop = False

            while (BUDGET > len(temp_best)):
                #percentage_done = round(np.floor((B0-BUDGET)/B0*1000)/10,2)
                #print(f"{percentage_done}% done")
                diversity = fn.get_last_diversity(temp_best)
                homogenity = 1 - diversity

                #crossover_ratio =  0.1#0.01 + ((diversity)/10)
                #print(f'crossover_ratio = {crossover_ratio}, N={N}')

                temp_pop, cost = fn.populate_from_df(temp_best,N,metric_function_list,mutation_function_list,
                                            generation+1,include_seeds=True,fitness='Metric 1',crossover=crossover, crossover_type=0,
                                            crossover_ratio=crossover_ratio, fitness_proportional = True)
                cost_history.append(cost)
                BUDGET -= cost
                temp_pop.reset_index(drop=True,inplace=True)


                #if homogenity>0:
                #    temp_best = fn.get_percent_best(temp_pop, metrics,0.05+((homogenity**2)/2),minimize=False)
                #else:
                temp_best = fn.get_percent_best(temp_pop, metrics,percent_best_fraction,minimize=False)
                gen_history = pd.concat([gen_history,temp_best])
                
                
                
                if include_initial_pop:
                    initial_sample_size = int(np.ceil(  len(temp_best) *  initial_sample_fraction   ))
                    #initial_sample_size = 3 * int(np.ceil( len(temp_best)* homogenity**2  ))
                    #if initial_sample_size>len(initial_pop):
                    #    initial_sample_size = len(initial_pop)

                    temp_best = pd.concat([temp_best, initial_pop.sample(initial_sample_size)])

            
                #if homogenity==0:
                N = int(np.ceil(BUDGET * next_generation_fraction))
                    #print(f"homo {homogenity} , N {N}")
                #else:
                #    N = int(np.ceil(BUDGET * 0.1  * homogenity))
                    #print(f"no homo N {N}")
                #print(f"------- {N} ------- BUDGET / GENERATIONS {BUDGET / GENERATIONS} * homogenity {homogenity} + 0.1  = {(BUDGET / GENERATIONS) * homogenity + 0.1}")
                if N<=int(np.ceil(minimal_next_generation_fraction*B0)):
                    N = int(np.ceil(minimal_next_generation_fraction*B0))

                generation += 1
                if generation>=max_gens:
                    break

            print(f"Done in {time.time()-t0}")
            if SAVE:
                gen_history.to_csv("out_exp/"+file_name+".csv")

###
latest_gen = gen_history[gen_history['Generation']==gen_history['Generation'].max()]
latest_gen.head()