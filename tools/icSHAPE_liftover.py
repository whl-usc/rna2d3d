#!/usr/bin/env python3

################################################################################
"""
Contact:    wlee9829@gmail.com 
Date:       2025_05_09 
Python:     python3.11
Script:     icSHAPE_liftover.py

Library strategy: icSHAPE

After demultiplexing, reads were stripped adaptors (PMID: 25411354) and mapped 
to the Homo Sapiens transcriptome - GRCh37.74 (hg19) using Bowtie2. All samples 
have a 13nt barcode on the 5'end which needs to be trimmed before mapping to 
the transcriptome. icSHAPE reverse transcriptase (RT) stops were isolated as 
the first nucleotide after the 5'end barcode. RT stops were used to calculated 
icSHAPE reactivity scores

Genome_build: Homo sapiens - GRCh37.74 (hg19)

Supplementary_files_format_and_content: 2 Tables (one each for In Vitro and In 
Vivo modified samples) of icSHAPE reactivities across human transcripts 
measured. Format: Transcript ID, Transcript Length (in nucleotides), RPKM 
value, icSHAPE reactivity scores per base for the length of each transcript (
NULL values denote bases not measured by the icSHAPE experiment)

Converts from TXT to bedgraph, then use the liftOver env to convert to hg38 
using the chain file.
"""
################################################################################
import csv
import gzip

# Function to parse GTF file and extract exon regions
def parse_gtf(gtf_file):
    exons = {}
    
    with gzip.open(gtf_file, "rt") as f:
        for line in f:
            fields = line.strip().split("\t")
            if fields[2] == "exon":  # Only process exon entries
                chrom = f"chr{fields[0]}"
                start = int(fields[3]) - 1  
                # GTF is 1-based, BEDGRAPH is 0-based
                end = int(fields[4])
                transcript_id = None

                # Extract transcript_id from attributes field
                for attribute in fields[8].split(";"):
                    if "transcript_id" in attribute:
                        transcript_id = attribute.split(" ")[2].strip('"')
                        break

                if transcript_id:
                    if chrom not in exons:
                        exons[chrom] = {}
                    if transcript_id not in exons[chrom]:
                        exons[chrom][transcript_id] = []
                    exons[chrom][transcript_id].append((start, end))
    return exons

# Function to read and clean reactivity data
def read_reactivity_data(input_file):
    reactivities = {}
    with open(input_file, "r") as infile:
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            transcript_id = row[0]
            values = row[3:]
            cleaned_values = [0 if v == "NULL" else float(v) for v in values]
            reactivities[transcript_id] = cleaned_values
    return reactivities

# Function to write BEDGRAPH file
def write_bedgraph(output_file, reactivities, exons):
    with open(output_file, "w") as outfile:
        writer = csv.writer(outfile, delimiter="\t")

        for chrom in exons:
            for transcript_id in exons[chrom]:
                if transcript_id not in reactivities:
                    continue  # Skip if no data for this transcript
                values = reactivities[transcript_id]
                index = 0
                for start, end in exons[chrom][transcript_id]:
                    for pos in range(start, end):
                        if index < len(values):
                            writer.writerow([chrom, pos, pos + 1, values[index]])
                            index += 1
                        else:
                            break

# Define input and output files
gtf_file = "Homo_sapiens.GRCh37.74.gtf.gz"
# Parse the GTF file to get exon ranges for each transcript
exons = parse_gtf(gtf_file)

input_file = "GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt"
output_file = "HS_293T_icSHAPE_InVivo_BaseReactivities.bedgraph"
reactivity_values = read_reactivity_data(input_file)
write_bedgraph(output_file, reactivity_values, exons)
print(f"BEDGRAPH file written to {output_file}")

input_file = "GSE74353_HS_293T_icSHAPE_InVitro_BaseReactivities.txt"
output_file = "HS_293T_icSHAPE_InVitro_BaseReactivities.bedgraph"
reactivity_values = read_reactivity_data(input_file)
write_bedgraph(output_file, reactivity_values, exons)
print(f"BEDGRAPH file written to {output_file}")

