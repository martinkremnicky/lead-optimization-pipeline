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
import seaborn as sns
import matplotlib.pyplot as plt
import os
import os


#print(os.getcwd())
# Changing the current working directory
if os.getcwd() == 'C:\\Users\\marti\\Desktop\\__SKOLA_VU\\__THESIS\\repo\\lead-optimization-pipeline':
    os.chdir(os.getcwd()+"\\pipeline")
#print(os.getcwd())

SEED_SIZE_LIST = [10, 100, 1000, 10000, 100000]
SEED_SIZE = SEED_SIZE_LIST[2]
SEED_LIST = [0,1,2]
SEED = SEED_LIST[2]
BUDGET_LIST = [100, 1000, 10000, 100000]
BUDGET = BUDGET_LIST[2]
GENERATIONS = 50




SAVE = True


celecoxib = 'O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'
troglitazone = 'O=C1NC(=O)SC1Cc4ccc(OCC3(Oc2c(c(c(O)c(c2CC3)C)C)C)C)cc4'
task_f = [partial(met.compound_similarity,target_smile = celecoxib),
         partial(met.compound_similarity,target_smile = troglitazone),
         partial(standard_benchmarks.median_camphor_menthol().objective.score_list),
         partial(standard_benchmarks.median_tadalafil_sildenafil().objective.score_list),
         partial(standard_benchmarks.hard_osimertinib().objective.score_list),
         partial(standard_benchmarks.zaleplon_with_other_formula().objective.score_list)]




task_names = ['celecoxib','troglitazone','med1','med2','osimertinib','zaleplon']

task_dict = dict(zip(task_names, task_f))

task_names_run = ['med2']#,'med2'] #BUG BEWARE, 1st function works, any other don't - use one function at a time

SEED_SIZE_LIST_run = [100, 1000, 10000, 100000]
BUDGET_LIST_run = [100, 1000, 10000]#,100000]


#f = standard_benchmarks.zaleplon_with_other_formula().objective.score_list()

metric_function_list = [
    #partial(standard_benchmarks.zaleplon_with_other_formula().objective.score_list)
    #partial(standard_benchmarks.median_tadalafil_sildenafil().objective.score_list)
#    partial(met.compound_similarity,target_smile = troglitazone)
]
metrics = ['Metric 1']


mutation_function_list = [
    partial(mut.replacement,n=1 ,ring_aware= False)
    ,partial(mut.addition,fragment_size=1,rings=False, branches=False)
    ,partial(mut.deletion,n=1)
]


def dosmn():
    for task_i in range(len(task_names_run)):
        task_start = time.time()
        curr_task_name = task_names_run[task_i]
        curr_task_f = task_dict[curr_task_name]
        print(f"Doing {curr_task_name}")
        metric_function_list = [curr_task_f]
        print(metric_function_list)
        current_dateTime = datetime.now()
        time_format = f"{current_dateTime.date()}_{str(current_dateTime.hour).zfill(2)}-{str(current_dateTime.minute).zfill(2)}-{str(current_dateTime.second).zfill(2)}"
        chem_lookup = {}

        try:
            pbar.close()
            print("Ignore the previous status bar. It is from the previous run.")
        except:
            pass


        for SEED in SEED_LIST[:]:#SEED_LIST[1:]:
            for SEED_SIZE in SEED_SIZE_LIST_run:

                recalculate = True

                t0 = time.time()
                print("reading and evaluating data...")
                if recalculate:
                    SEED_PATH = f"../data/seed_{SEED}/rand_{SEED_SIZE}.tsv"
                    seed_df = pd.read_table(SEED_PATH) 
                    initial_pop, _, chem_lookup = fn.seeds_to_pop(seed_df,metric_function_list)
                    print('--',metric_function_list)
                else:
                    INITIAL_PATH = f"../data/evaluated_datasets/{curr_task_name}_{SEED}_{SEED_SIZE}.csv"
                    initial_pop = pd.read_csv(INITIAL_PATH)
                    chem_lookup = fn.dataframe_to_dict(initial_pop)
                print(f"Finished reading in {time.time()-t0}")


                #chem_lookup.update(initial_pop.set_index('SMILES molecule')['Metric 1'].to_dict()) #canon and set
                for BUDGET in BUDGET_LIST_run:
                    setup_name = f"seed_{SEED}_rand_{SEED_SIZE}"
                    file_name = f"{time_format}_{setup_name}_budget_{BUDGET}"
                    B0 = copy(BUDGET)
                    generation = 0
                    t0 = time.time()
                    
                    fit = 'Fitness'


                    first_generation_fraction = 0.05 #hyperparam
                    N = int(BUDGET * first_generation_fraction) #1st generation size

                    if len(initial_pop)>=N:
                        initial_best = initial_pop.sort_values(fit,ascending=False).head(N)
                    else:
                        initial_best = initial_pop
                    gen_history = pd.DataFrame(initial_pop)
                    gen_best_history = pd.DataFrame(initial_best)
                    temp_best = initial_best.copy(deep=True)
                    cost_history = []

                    #rest of hyperparams
                    crossover = True
                    crossover_ratio = 0.2

                    decreasing_size = True #if True, (leftover) BUDGET * next_generation_fraction = new N, else B0 * nfg = N
                    next_generation_fraction = 0.05 #size of next generation

                    minimal_next_generation_fraction = 0.0025
                    max_gens = 500
                    stagnation_gens = 75

                    canonical_best = False
                    percent_best_fraction = 0.1

                    include_initial_pop = False
                    initial_sample_fraction = 0.5 #initial_sample_size = len(temp_best) *  initial_sample_fraction ))

                    diversity_decay = True
                    selection_exponenet = 1


                    pbar = tqdm(desc=f"Working with budget of {B0}",
                                total=B0)
                    while (BUDGET > len(temp_best)):
                        #diversity = fn.get_last_diversity(temp_best)
                        #homogenity = 1 - diversity
                        temp_pop, cost, chem_lookup = fn.populate_from_df(temp_best,N,metric_function_list,mutation_function_list,
                                                    generation+1,include_seeds=True,fitness=fit,crossover=crossover, crossover_type=0,
                                                    crossover_ratio=crossover_ratio, proportional_sampling = True, 
                                                    lookup_dict=copy(chem_lookup), randomize_seeds=True,
                                                    avg_dist=diversity_decay,
                                                    exp=selection_exponenet,
                                                    use_fitness = True)
                        cost_history.append(cost)
                        BUDGET -= cost
                        temp_pop.reset_index(drop=True,inplace=True)
                        gen_history = pd.concat([gen_history,temp_pop])
                        gen_history.reset_index(drop=True,inplace=True)

                        if canonical_best:
                            temp_best = fn.get_canonical_percent_best(temp_pop,[fit],percent_best_fraction,minimize=False)
                        else:
                            temp_best = fn.get_percent_best(temp_pop,[fit],percent_best_fraction,minimize=False)    
                        gen_best_history = pd.concat([gen_best_history,temp_best])
                        gen_best_history.reset_index(drop=True,inplace=True)
                        
                        
                        
                        if include_initial_pop:
                            initial_sample_size = int(np.ceil(  len(temp_best) *  initial_sample_fraction   ))
                            temp_best = pd.concat([temp_best, initial_pop.sample(initial_sample_size)])
                        
                        if decreasing_size:
                            N = int(np.ceil(BUDGET * next_generation_fraction))
                        else:
                            N = int(np.ceil(B0 * next_generation_fraction))
                        if N<=int(np.ceil(minimal_next_generation_fraction*B0)):
                            N =int(np.ceil(minimal_next_generation_fraction*B0))
                        pbar.update(cost)
                        generation += 1

                        if generation>=max_gens:
                            break

                        
                        #current_gen = temp_pop['Generation'].max()


                        if generation >= stagnation_gens:
                            current_metric_max = temp_pop[temp_pop['Generation'] == generation]['Metric 1'].max()
                            metric_max_stagnation_gens_ago = gen_best_history[gen_best_history['Generation'] == (generation - stagnation_gens)]['Metric 1'].max()
                            print(f"COMAPRE {current_metric_max} VS OLD {metric_max_stagnation_gens_ago} - {generation- stagnation_gens}")
                            if current_metric_max <= metric_max_stagnation_gens_ago:
                                print(f"#Stagnation for more than {stagnation_gens} generations")
                                break



                    pbar.update(BUDGET)
                    pbar.close()
                    print(f"Done {generation} in {time.time()-t0}")

                    exp_path = "exp_24"


                    FOLDER_P = f"{exp_path}/best_{curr_task_name}"
                    if SAVE:
                        if not os.path.exists(FOLDER_P):
                            os.makedirs(FOLDER_P)
                        gen_best_history.to_csv(f"{FOLDER_P}/{file_name}.csv")
                    FOLDER_P = f"{exp_path}/all_{curr_task_name}"
                    if False:
                        if not os.path.exists(FOLDER_P):
                            os.makedirs(FOLDER_P)
                        gen_history.to_csv(f"{FOLDER_P}/{file_name}.csv")
        print(f"!!! ---  Done with {curr_task_name} in {time.time() - task_start}")
        ###
        latest_gen = gen_best_history[gen_best_history['Generation']==gen_best_history['Generation'].max()]
        pd.set_option('display.max_columns', None)
        print(latest_gen.head())




if __name__ == '__main__':
    import cProfile, pstats
    profiler = cProfile.Profile()
    profiler.enable()
    dosmn()
    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats('tottime')
    stats.print_stats(50)
    print("===================================================\n"*3)
    stats = pstats.Stats(profiler).sort_stats('cumtime')
    stats.print_stats(50)