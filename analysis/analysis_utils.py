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
            
            # Get the median of "Metric" where "Generation" is 10
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
def plot_medians(medians_dict):
    # Loop over the sizes and plot the medians
    for size, medians in medians_dict.items():
        budgets = sorted(medians.keys())
        median_values = [medians[budget] for budget in budgets]
        
        # Plot on a log scale
        plt.plot((budgets), median_values, label=f"IDS: {size}")
    
    plt.xlabel("Budget (log scale)")
    plt.ylabel("Median Score")
    plt.title("Median Best Guacamol 'Median Molecule 1' Score of 3 runs,\nby Budget and Initial Dataset Size (IDS)")
    plt.xscale('log')
    plt.legend()
    plt.show()

def output_coordinates(medians_dict):
    output_strings = []
    
    # Loop over each size and its budgets and medians
    for size, budgets_medians in medians_dict.items():
        # Get the budgets and medians in ascending order of budget
        budgets_medians_sorted = sorted(budgets_medians.items())
        
        # Convert the budgets and medians to strings and format them as coordinates
        coordinates = " ".join(f"({budget}, {round(median,3)})" for budget, median in budgets_medians_sorted)
        # Format the size and coordinates as a string and add it to the list
        output_strings.append("                \\addplot coordinates {"+f"{coordinates}"+"}; "+f"%dataset size = {size}")
    
    beginning = '    \\begin{subfigure}[b]{0.45\\textwidth}\n        \\begin{tikzpicture}\n            \\begin{axis}[xmode=log, ymin=0.45, ymax=1, width=\linewidth,\n            cycle list={cyan!100!lime, cyan!66!lime, cyan!33!lime, cyan!0!lime}]'
    end = '            \end{axis}\n        \\end{tikzpicture}\n        %\\caption{Plot 1}\n    \\end{subfigure}'
    # Join all the strings with line breaks and return the result
    return beginning+"\n"+"\n".join(output_strings)+'\n'+end




# Process the files and plot the medians
medians_dict = process_files("./pipeline/out_exp/")
coordinates_string = output_coordinates(medians_dict)
print(coordinates_string)
plot_medians(medians_dict)

