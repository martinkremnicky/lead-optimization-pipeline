import pandas as pd
import os

def add_columns_from_filename(file_path):
    # Extract the filename from the file_path
    base_name = os.path.basename(file_path)
    
    # Strip the .csv from the filename
    file_name = os.path.splitext(base_name)[0]
    
    # Split the filename into its components
    filename_parts = file_name.split("_")
    
    # Extract seed, size, and budget from the filename_parts
    seed = filename_parts[filename_parts.index('seed')+1]
    size = filename_parts[filename_parts.index('rand')+1]
    budget = filename_parts[filename_parts.index('budget')+1]
    
    # Read the csv file into a pandas DataFrame
    df = pd.read_csv(file_path, index_col=0)
    
    # Add new columns to the dataframe
    df['seed'] = int(seed)
    df['size'] = int(size)
    df['budget'] = int(budget)
    
    # Return the modified DataFrame
    return df

# Test the function
#df = add_columns_from_filename("./pipeline/out/2023-06-12_17-05-06_seed_0_rand_1000_budget_1000.csv")
#print(df.head())


import pandas as pd
import os
import matplotlib.pyplot as plt
import glob

def process_files(folder_path):
    # Collect all .csv files in the given folder
    files = glob.glob(os.path.join(folder_path, "*.csv"))
    
    # Create a dictionary to store data from each file
    data = {size: {} for size in set([get_size(file) for file in files])}
    
    # Process each file
    for file_path in files:
        df = add_columns_from_filename(file_path)
        median = df[df["Generation"] == 10]["Metric 1"].median()
        
        # Get seed, size and budget from filename
        size = get_size(file_path)
        budget = get_budget(file_path)
        
        # Store the median
        if budget in data[size]:
            data[size][budget].append(median)
        else:
            data[size][budget] = [median]

    # For each size, average the medians and plot
    for size, budgets in data.items():
        budgets = {k: sum(v)/len(v) for k, v in budgets.items()}  # Average medians
        sorted_budgets = dict(sorted(budgets.items()))  # Sort by budget
        plt.plot(list(sorted_budgets.keys()), list(sorted_budgets.values()), label=f"Size: {size}")

    # Set plot labels and title
    plt.xlabel('Budget')
    plt.ylabel('Median of Metric')
    plt.title('Median of Metric for Generation 10 per Budget and Size')
    plt.legend()
    plt.show()

# Helper functions to get size and budget from file name
def get_size(file_path):
    base_name = os.path.basename(file_path)
    file_name = os.path.splitext(base_name)[0]
    filename_parts = file_name.split("_")
    return int(filename_parts[filename_parts.index('rand')+1])

def get_budget(file_path):
    base_name = os.path.basename(file_path)
    file_name = os.path.splitext(base_name)[0]
    filename_parts = file_name.split("_")
    return int(filename_parts[filename_parts.index('budget')+1])

# Call the function with the folder where your csv files are located
#process_files("./pipeline/out/")
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

def process_files(directory):
    # Initialize a dictionary to hold the medians
    medians_dict = {}
    
    # Loop over all files in the directory
    for file_name in os.listdir(directory):
        if file_name.endswith(".csv"):
            file_path = os.path.join(directory, file_name)
            
            # Extract the necessary data from the filename
            file_name = file_name[:-4] #.csv
            filename_parts = file_name.split("_")
            seed = int(filename_parts[filename_parts.index('seed')+1])
            size = int(filename_parts[filename_parts.index('rand')+1])
            budget = int(filename_parts[filename_parts.index('budget')+1])
            
            # Read the csv file into a pandas DataFrame
            df = pd.read_csv(file_path)
            
            # Get the median of "Metric" where "Generation" is 10
            median = df[df["Generation"] == df["Generation"].max()]["Metric 1"].median()
            
            # Add the median to the dictionary
            if size not in medians_dict:
                medians_dict[size] = {}
            medians_dict[size][budget] = median
            
    return medians_dict

import matplotlib as mpl
import numpy as np

# Set the default color cycle
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#d62728','#ff7f0e','#2ca02c','#9467bd']) 

def plot_medians(medians_dict):
    # Loop over the sizes and plot the medians
    for size, medians in medians_dict.items():
        budgets = sorted(medians.keys(), reverse=True)
        median_values = [medians[budget] for budget in budgets]
        
        # Plot on a log scale
        plt.plot((budgets), median_values, label=f"Initial dataset size {size}")
    
    plt.xlabel("Budget (log scale)")
    plt.ylabel("Median Metric")
    plt.title("Median Metric by Budget and Size")
    plt.xscale('log')
    plt.legend()
    plt.show()

# Process the files and plot the medians


import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def process_files(directory):
    # Initialize a dictionary to hold the medians
    medians_dict = defaultdict(lambda: defaultdict(list))
    
    # Loop over all files in the directory
    for file_name in os.listdir(directory):
        if file_name.endswith(".csv"):
            file_path = os.path.join(directory, file_name)
            file_name = file_name[:-4] #.csv
            
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
    plt.title("Median Best Guacamol Zaleplon Score of 3 runs,\nby Budget and Initial Dataset Size (IDS)")
    plt.xscale('log')
    plt.legend()
    plt.show()

# Process the files and plot the medians
medians_dict = process_files("./pipeline/out_zaleplon_all/")
plot_medians(medians_dict)

