"""
contact:    wlee9829@gmail.com
date:       2024_05_10
python:     python3.10
script:     ribomethscore.py

This script is an updated version of the ribomethscore_bed2bedgraph.py, 
used to calculate the RiboMethScore(ScoreC2).

"""
# Define version
__version__="1.1.0"

# Version notes
__update_notes__="""
1.1.0
    -   Added optional genelist argument to filter reads, faster calculations.
    -   Fixed function formatting for readability.

1.0.0
    -   Initial commit, no functional changes to original script. 
    -   Adjusted formatting for clarity.
"""

# Import Packages
from datetime import datetime
import argparse
import itertools
from itertools import combinations
import matplotlib.pyplot as plt 
from multiprocessing import Process, Lock, Manager
import numpy
import os
import random
import re
import subprocess
import sys
import textwrap
import time

###########################################################################
# 1. Define the functions.

def timenow():
    """
    Returns the current timestamp as a string.

    Returns:
        str: Current timestamp in format 'YYY-MM-DD HH:MM:SS'.
    """
    time = str(datetime.now())[:-7]

    return time

def readGenefile(genebed):
    """
    Reads the genes bedfile into a dictionary.

    Returns:
        dict: gene_start, gene_end, and gene_name.
    """
    print(str(datetime.now())[:-7], "Reading gene bed file...")
    gene_dict = {}    
    with open(genebed, 'r') as f:
        for line in f:
            align = line.rstrip("\n").split("\t")
            gene_start, gene_end, gene_name =\
                int(align[1]), int(align[2]), str(align[3])
            gene_dict[gene_name] = gene_end - gene_start + 1

    return gene_dict

def intersect(alignbed, genebed, output):
    """
    Reads the sorted alignments bedfile into a dictionary, uses bedtools
    intersect to determine alignment and gene bedfile overlaps.

    Returns:
        bedgraph: 
    """
    print(str(datetime.now())[:-7], "Reading alignment bed file, calculating"
        "overlap...")

    # Intersect the alignment and gene bedfiles.
    output_tmp = output + '_overlap.bed'
    intersect = ("bedtools intersect -a %s -b %s -wa -wb -sorted -f 0.5 > %s"\
        % (alignbed, genebed, output_tmp))
    os.system(intersect)

    # Process the overlap bed file
    bedgraph = {}; gene_cov = {}; 
    with open(output_tmp, 'r') as f:
        for line in f:
            align = line.rstrip('\n').split('\t')
            
            chr_name, read_start, read_end, \
            gene_start, gene_end, gene_name, gene_strand = \
            align[0], int(align[1]) + 1, int(align[2]), \
            int(align[7]), int(align[8]), align[9], align[11]
        
            if gene_name not in gene_cov:
                gene_cov[gene_name] = {}

            # 5'end
            if read_start >= gene_start and read_start <= gene_end:
                chr_pos = chr_name+'_'+str(read_start)
                if chr_pos not in bedgraph:
                    bedgraph[chr_pos] = 0
                bedgraph[chr_pos] += 1
                #count the coverage based on gene
                if gene_strand == "+": 
                    gene_pos = read_start - gene_start + 1
                if gene_strand == "-": 
                    gene_pos = gene_end - read_end + 1
                gene_pos = gene_name+'_'+str(gene_pos)
                if gene_pos not in gene_cov[gene_name]:
                    gene_cov[gene_name][gene_pos] = 0
                gene_cov[gene_name][gene_pos] += 1
                    
            # 3'end
            if read_end >= gene_start and read_end <= gene_end:
                chr_pos = chr_name+'_'+str(read_end)
                if chr_pos not in bedgraph:
                    bedgraph[chr_pos] = 0
                bedgraph[chr_pos] += 1
                #count the coverage based on gene
                if gene_strand == "+":
                    gene_pos = read_end - gene_start + 1
                if gene_strand == "-":
                    gene_pos = gene_end - read_end + 1
                gene_pos = gene_name+'_'+str(gene_pos)
                if gene_pos not in gene_cov[gene_name]:
                    gene_cov[gene_name][gene_pos] = 0
                gene_cov[gene_name][gene_pos] += 1

    # output bedgraph file
    with open(output+'_5_3end.bedgraph', 'w') as output_bedgraph:
        for chr_pos in sorted(bedgraph):
            chr_name, pos = chr_pos.split('_')
            start_pos = int(pos) - 1
            output_bedgraph.write(f"{chr_name}\t{start_pos}\t{pos}"
                f"\t{bedgraph[chr_pos]}\n")

    return gene_cov

def RiboMethScore2(l2,l1,x,r1,r2):
    """
    Calls the RiboMethScore2 (ScoreC2) for 2 positions left and right of the
    specified target position.

    l2: reads count at position l2 (2nd left flanking regions of i)
    l1: reads count at position l1 (1st left flanking regions of i)
    x:  reads count at target position
    R1: reads count at position r1 (1st right flanking regions of i)
    R2: reads count at position r2 (2nd right flanking regions of i)

    Returns:
        int: C2 score
    """    
    if all (l == 0 for l in [l2, l1]) and (r == 0 for r in [r1, r2]):
        ScoreC2 = 0
    else:
        ScoreC2 = 1-x/(0.5*((l2*0.9+l1*1.0)/(1.0+0.9) + 
            (r1*1.0+r2*0.9)/(1.0+0.9)))
    if ScoreC2 > 0:
        return ScoreC2
    else:
        return 0

def RiboMethScore1(l1,x,r1):
    """
    Calls the RiboMethScore6 (ScoreC1) for 1 position left and right of the 
    specified target position.

    L1: reads count at position l1 (Y left flanking regions of i)
    x:  reads count at target position
    R1: reads count at position ry (Y right flanking regions of i)
 
    Returns:
        int: C1 score
    """
    if l1==0 and r1==0:
        ScoreC1 = 0
    else:
        ScoreC1 = 1-x/(0.5*( (l1*1.0)/(1.0) + (r1*1.0)/(1.0) ))
    if ScoreC1 > 0:
        return ScoreC1
    else:
        return 0

def calculateRMS2(gene_dict, gene_cov, output, gene_list=False):
    """
    Calculates the RiboMethScore2.
    """
    print(str(datetime.now())[:-7], "Calculate RiboMethScore2 ...")
    if gene_list:
        with open('genelist.txt', 'r') as f:
            genes_list = [line.rstrip("\n") for line in f]

    with open(output+'_Methscore.txt', 'w') as output_ribomethscore: 
        for gene in gene_cov:
            if not gene_list or gene in genes_list:
                for i in range(3, gene_dict[gene] - 1, 1):
                    l2 = 0; l1 = 0; x = 0; r1 = 0; r2 = 0
                    if gene+'_'+str(i-2) in gene_cov[gene]:
                        l2 = gene_cov[gene][gene+'_'+str(i-2)]
                    if gene+'_'+str(i-1) in gene_cov[gene]:
                        l1 = gene_cov[gene][gene+'_'+str(i-1)]
                    if gene+'_='+str(i) in gene_cov[gene]: 
                        x = gene_cov[gene][gene+'_'+str(i)]
                    if gene+'_'+str(i+1) in gene_cov[gene]: 
                        r1 = gene_cov[gene][gene+'_'+str(i+1)]
                    if gene+'_'+str(i+2) in gene_cov[gene]: 
                        r2 = gene_cov[gene][gene+'_'+str(i+2)]
                    string = gene+'_'+str(i)+'\t'+\
                        str(RiboMethScore2(l2,l1,x,r1,r2))+\
                        '\t'+str((l1+x+r1)/3)+'\n'
                    output_ribomethscore.write(string)

############################################################################

# 2. Main, define the accepted arguments.
def main():
    parser = argparse.ArgumentParser(
        prog="ribomethscore.py",
        description=textwrap.dedent("""\
###########################################################################
Pass
###########################################################################
"""),
usage='''

- For calculating the RiboMeth score and writing out to bedgraph file.
\npython3 %(prog)s
''')    
    parser.add_argument("alignbed", help="Bed file from aligned bam file.")
    parser.add_argument("genebed", help="Gene bed information.")
    parser.add_argument("--genelist", "-g", help="Gene list to calculate RMS2"
        " for.", default=[])
    parser.add_argument("output", help="Output filename.")

    args = parser.parse_args()

    # Check required arguments
    genebed = args.genebed
    alignbed = args.alignbed
    genelist = args.genelist
    output = args.output

    gene_dict = readGenefile(genebed)
    gene_cov = intersect(alignbed, genebed, output)
    calculateRMS2(gene_dict, gene_cov, output, genelist)

if __name__ == "__main__":
    main()
sys.exit()