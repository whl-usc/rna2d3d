"""
contact:    wlee9829@gmail.com
date:       2024_07_03
python:     python3.10
script:     plot_readcount.py

This Python script plots gene counts in user-defined styles. 
Graph can be log(2) transformed to find normal distribution.
"""
# Define version
__version__ = "1.0.0"

# Version notes
__update_notes__ = """
1.0.0   
    -   Initial commit, set up functions.
"""

# Import Packages
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

###########################################################################
def readfile(file_path, cutoff=None, exclude_chromosomes=None):
    """
    Opens coverage text file and parses information.
    """
    column_names = ['chrom', 'start', 'end', 'gene_name', 
        'score', 'strand', 'count']
    
    df = pd.read_csv(file_path, delimiter='\t', 
        header=None, names=column_names)

    extracted_df = df[['chrom', 'gene_name', 'count']].copy()

    if cutoff is not None:
        extracted_df = extracted_df[extracted_df['count'] > cutoff]

    if exclude_chromosomes is not None:
        filtered_df = extracted_df[~extracted_df['chrom'].isin(exclude_chromosomes)]
    else:
        filtered_df = extracted_df

    return filtered_df

def plot_gene_counts(dataframes, output_file, log_transform=False):
    """
    Plots histogram of gene counts from the dataframe, optionally
    log2-transformed.
    """
    plt.figure(figsize=(12, 6))  

    for dataframe, filename in dataframes:
        if log_transform:
            dataframe['log2_count'] = np.log2(dataframe['count'] + 1)
            count_data = dataframe['log2_count']
            xlabel = 'Log2 Transformed Expression Counts'
            title = 'Overlay of Log2 Transformed Gene Counts'
            bin_count = int(np.ceil(np.sqrt(len(dataframe))))
        else:
            count_data = dataframe['count']
            xlabel = 'Raw Expression Counts'
            title = 'Overlay of Raw Gene Counts'
            plt.xlim(0, 20000)
            bin_count = 2000

        plt.hist(count_data, bins=bin_count, edgecolor='black', 
            alpha=0.5, label=filename)

    plt.xlabel(xlabel)
    plt.ylabel('Frequency (Number of Genes)')
    plt.title(title)
    plt.legend(loc='upper right')
    plt.grid(False)
    plt.savefig(f"{output_file}.png", dpi=400)
    plt.close()

def parse_args():
    """
    Main function to set up argument parsing, handles input arguments, calls
    relevant functions for processing and plotting read count data.
    """
    parser = argparse.ArgumentParser(
        prog="plot_readcount.py",        
        description="Plot gene counts from tab-separated files.")
    
    parser.add_argument("input_files", type=str, nargs='+',
        help="Paths to the input files containing gene count data.")

    parser.add_argument("output_file", type=str,
        help="Path to save the output plot as a PNG file.")

    parser.add_argument("--cutoff", type=int, default=None,
        help="Cutoff value for filtering gene counts.")

    parser.add_argument("--log2", action='store_true',
        help="Transform gene counts using log2.")

    return parser.parse_args()

def main(args):
    input_files = args.input_files
    output_file = args.output_file
    cutoff = args.cutoff
    log2_transform = args.log2

    exclude_chromosomes = None
    if not log2_transform:
        exclude_chromosomes = ["mm12S", "mm16S", "mm45S", "mm4.5S", 
                               "mm5S", "mmBC1", "mmsnRNA", "RN7SK", 
                               "RN7SL", "RNU7", "RNY"]

    dataframes = []
    for file in input_files:
        data = readfile(file, cutoff, exclude_chromosomes)
        filename = file.split(".bedtools_coverage.txt")[0]
        dataframes.append((data, filename))
        print(filename)
        
    plot_gene_counts(dataframes, output_file, log2_transform)

if __name__ == "__main__":
    args = parse_args()
    main(args)
    sys.exit()