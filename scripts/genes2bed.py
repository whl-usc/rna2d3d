'''
contact:    wilsonhl@usc.edu
date:       2022_09_14
python:     python3.10
script:     genes2bed.py
    
This script is used to generate a bed file for CRSSANT analysis.
Contains the start, end of genes over a specific number of
reads as determined by the CountCdsUtr_exclude.py script. 
'''

# Import packages
import sys, argparse, os, re, itertools, math, time
import pandas as pd
from datetime import datetime

# Usage instructions
if len(sys.argv) < 5:
    print("Usage:           python genes2bed.py *_ReadCount.txt CdsUtr.bed min_reads outname")
    print("*_ReadCount.txt: Output from CountCdsUtr_exclude.py, shows read counts")
    print("CdsUtr.bed:      Bed file of CDS, 5UTR, 3UTR regions and strand")
    print("min_reads:       Minimum amount of reads before the gene is included")
    print("outname:         Output file prefix")
    sys.exit()

in_file = open(sys.argv[1],'r')
bed_file = open(sys.argv[2], 'r')
min_reads = int(sys.argv[3])
outname = str(sys.argv[4])

if min_reads < 0:
    print("Read count values must be >0")    
    sys.exit()

####################################################################################################

# Process the text file into a dataframe 
df = pd.read_csv(in_file, sep="\t")
df_minread = df.loc[((df['gene_Reads'] >= min_reads) & (df['Biotype'] =='protein_coding'))]
genes_list = df_minread['Gene'].values.tolist(); num_reads=len(genes_list)
print(str(datetime.now())[:-7] + " Searching *_ReadCount.txt based on min_reads: " + str(min_reads))
# Print the list of genes to a text file for further analysis.
# print(df_genes, file=open(outname+'.txt','a'))

# Open CdsUtr.bed file for parsing.
print("                 Searching " + str(sys.argv[2]) + " for: " + str(num_reads) + " genes")
bed = pd.read_csv(bed_file, sep="\t")
bed.columns=["chromosome","start","stop","strand","gene","type","class"]

# Collapse rows for each gene and write to dictionary.
print(str(datetime.now())[:-7] + " Compiling genes to bed file...")
rows = []; pc=bed.reset_index()
for m,n in enumerate(genes_list):
    a = pc.loc[(pc['gene'] == n)]
    chromosome = a.iloc[0,1]
    start = str(a.iloc[0,2])
    stop = str(a.iloc[-1,3])
    gene = n
    coverage = str(1000)
    strand = a.iloc[0,4]
    row = {'chromosome':chromosome, 'start':start, 'stop':stop, 'gene:':gene, 'coverage':coverage, 'strand':strand} 
    #bed_out = bed_out.append(pd.DataFrame([rows]), ignore_index=True)
    rows.append(row)
    if (m+1) % 1000 == 0:
        print("                   " + " Processed " + str(m+1) + " genes.")

bed = pd.DataFrame(rows); outfile = outname+".bed"
bed.to_csv(outfile, sep='\t',index=False, header=None)
print(str(datetime.now())[:-7] + " Completed; see " + outfile)
