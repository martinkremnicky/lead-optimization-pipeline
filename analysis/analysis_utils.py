import pandas as pd
import os


import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import matplotlib as mpl

def process_files(directory):
    # Initialize a dictionary to hold the medians
    medians_dict = defaultdict(lambda: defaultdict(list))
    
    # Loop over all files in the directory
    for file_name in os.listdir(directory):
        if file_name.endswith(".csv"):
            file_path = os.path.join(directory, file_name)
            file_name = file_name[:-4] #delete '.csv'

            
            # Extract the necessary data from the filename
            filename_parts = file_name.split("_")
            seed = int(filename_parts[filename_parts.index('seed')+1])
            size = int(filename_parts[filename_parts.index('rand')+1])
            budget = int(filename_parts[filename_parts.index('budget')+1])
            
            # Read the csv file into a pandas DataFrame
            df = pd.read_csv(file_path)
            
            # Get the median of "Metric 1" where "Generation" is max
            median = df[df["Generation"] == df["Generation"].max()]["Metric 1"].median()
            
            # Add the median to the dictionary
            medians_dict[size][budget].append(median)
    
    # Get the overall median for each combination of size and budget
    for size, budgets in medians_dict.items():
        for budget, medians in budgets.items():
            medians_dict[size][budget] = np.median(medians)
            
    return medians_dict
['#9467bd','#2ca02c','#ff7f0e','#d62728']
['#d62728','#ff7f0e','#2ca02c','#9467bd']
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#9467bd','#2ca02c','#ff7f0e','#d62728']) 
def plot_medians(medians_dict, task_name):
    # Loop over the sizes and plot the medians
    for size, medians in medians_dict.items():
        budgets = sorted(medians.keys())
        median_values = [medians[budget] for budget in budgets]
        
        # Plot on a log scale
        plt.plot((budgets), median_values, label=f"IDS: {size}")
    
    plt.xlabel("Budget (log scale)")
    plt.ylabel("Median Score")
    plt.title(f"Median Best Guacamol {task_name} Score of 3 runs,\nby Budget and Initial Dataset Size (IDS)")
    plt.xscale('log')
    plt.legend()
    plt.show()

def output_coordinates(medians_dict,task_name:str):
    output_strings = []
    max_median = 0
    # Loop over each size and its budgets and medians
    print(task_name)
    t_n_s = task_name.split("_")
    for size, budgets_medians in medians_dict.items():
        # Get the budgets and medians in ascending order of budget
        budgets_medians_sorted = sorted(budgets_medians.items())
        #print(budgets_medians_sorted)
        for budget, median in budgets_medians_sorted:
            if median> max_median:
                max_median = round(median*1.1,3) #JUST FOR THE PLOT IN LaTeX
        # Convert the budgets and medians to strings and format them as coordinates
        coordinates = " ".join(f"({budget}, {round(median,3)})" for budget, median in budgets_medians_sorted)
        # Format the size and coordinates as a string and add it to the list

        page_frac = 0.45
        for budget, median in budgets_medians_sorted:
            print(f"{size}")
        for budget, median in budgets_medians_sorted:
            print(f"{budget}")
        for budget, median in budgets_medians_sorted:
            print(f"{median}")
 
        color_settings = 'cycle list={cyan!100!lime, cyan!66!lime, cyan!33!lime, cyan!0!lime}'


        output_strings.append("                \\addplot coordinates {"+f"{coordinates}"+"}; "+f"%dataset size = {size}")
    
    beginning = '    \\begin{subfigure}[b]{'+str(page_frac)+'\\textwidth}\n        \\begin{tikzpicture}\n            \\begin{axis}[ytick=\\empty, xmode=log, ymin=0.00, ymax='+f"{max_median}"+', width=\linewidth,\n            '+color_settings+']'
    end = '            \end{axis}\n        \\end{tikzpicture}\n        %\\caption{'+task_name+'}\n    \\end{subfigure}'
    # Join all the strings with line breaks and return the result
    return beginning+"\n"+"\n".join(output_strings)+'\n'+end



Z = 0

#task_names = ['celecoxib','troglitazone','med1','med2','osimertinib','zaleplon']
task_names = ['troglitazone','med2']
#task_names = ['med2']
for i in range(1,1+14):
    pass
exp = 26

#for i in range(1,1+26):
#    exp = i
if True:
    for task_n in task_names[:]:
        try:
            # Process the files and plot the medians
            medians_dict = process_files(f"./pipeline/exp_{exp}/best_{task_n}/")
            #medians_dict = process_files(f"./pipeline/old_results/out_{task_n}/")
            coordinates_string = output_coordinates(medians_dict,f"exp_{exp}_{task_n}")
            print()
            print(f"%exp_{exp}_{task_n} below")
            print(coordinates_string)
            print("\n"*5)
            plot_medians(medians_dict, f"exp_{exp}_{task_n}")
        except:
            pass

