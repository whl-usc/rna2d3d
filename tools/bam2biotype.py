#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_03_04
Python:     python3.10
Script:     bam2biotypes.py
"""

################################################################################
# Define version
__version__ = "1.0.0"

# Version notes
__update_notes__ = """
1.0.0
    - Initial commit, using gffutils for GTF parsing.
"""

################################################################################
# Import packages
import argparse
import os
import pysam
import pandas as pd
from collections import defaultdict

################################################################################
# Define sub-functions for processing

def parse_gene(gtf_file):
    """
    Parses a GTF file into a DataFrame, extracting gene_id and gene_biotype.
    """
    nc_to_chr = {
        "NC_000001.11": "chr1", "NC_000002.12": "chr2", "NC_000003.12": "chr3",
        "NC_000004.12": "chr4", "NC_000005.10": "chr5", "NC_000006.12": "chr6",
        "NC_000007.14": "chr7", "NC_000008.11": "chr8", "NC_000009.12": "chr9",
        "NC_000010.11": "chr10", "NC_000011.10": "chr11", 
        "NC_000012.12": "chr12", "NC_000013.11": "chr13", 
        "NC_000014.9": "chr14", "NC_000015.10": "chr15",
        "NC_000016.10": "chr16", "NC_000017.11": "chr17", 
        "NC_000018.10": "chr18", "NC_000019.10": "chr19", 
        "NC_000020.11": "chr20", "NC_000021.9": "chr21", 
        "NC_000022.11": "chr22", "NC_000023.11": "chrX", "NC_000024.10": "chrY"
    }

    col_names = ([
        "chrom", "source", "feature", 
        "start", "end", "score", 
        "strand", "frame", "attributes"])

    df = pd.read_csv(
        gtf_file, sep="\t", comment="#", names=col_names, dtype=str)

    df = df[df["feature"] == "gene"]
    df["chrom"] = df["chrom"].replace(nc_to_chr)

    df["gene_id"] = df["attributes"].str.extract(r'gene_id "([^"]+)"')
    df["gene_biotype"] = df["attributes"].str.extract(r'gene_biotype "([^"]+)"')
    df = df[["chrom", "start", "end", "strand", "gene_id", "gene_biotype"]]

    # Define the custom genes, each with its own chromosome name
    custom_genes = pd.DataFrame({
        "chrom": [
            "hs12S", "hs16S", "hs5S", "hs45S", "RN7SK", "RN7SL", "hssnRNA",
            "RNU7", "RNY", "U3", "U8", "U13", "U14AB", "U17"
        ],
        "start": ["1"] * 14,
        "end": [
            "954", "1559", "121", "13357", "331", "288", "2058", "63",
            "694", "217", "136", "120", "283", "207"
        ],
        "strand": ["+"] * 14,
        "gene_id": [
            "hs12S", "hs16S", "hs5S", "hs45S", "RN7SK", "RN7SL", "hssnRNA",
            "RNU7", "RNY", "U3", "U8", "U13", "U14AB", "U17"
        ],
        "gene_biotype": [
            "rRNA", "rRNA", "rRNA", "rRNA", "snRNA", "snRNA", "snRNA",
            "snRNA", "Y_RNA", "snoRNA", "snoRNA", "snoRNA", "snoRNA",
            "snoRNA"
        ]
    })

    # Append custom genes to the parsed GTF data
    df = pd.concat([df, custom_genes], ignore_index=True)
    return df


def count_biotype(df, bamfile):
    """
    Counts reads assigned to each gene biotype using BAM file.
    Returns a dictionary with biotype counts for exon and intron reads.
    """
    biotype_counts = defaultdict(int)

    # Load BAM file
    bam = pysam.AlignmentFile(bamfile, "rb")
    bam_chroms = set(bam.references)
    
    # Iterate over genes
    for _, row in df.iterrows():
        chrom, start, end, gene_biotype = row["chrom"], int(row["start"]), int(row["end"]), row["gene_biotype"]

        if chrom not in bam_chroms:
           continue

        # Count reads overlapping this gene
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped:
                continue
            biotype_counts[gene_biotype] += 1

    bam.close()
    return biotype_counts


def main():
    """Main function for processing BAM file and summarizing biotype counts."""
    parser = argparse.ArgumentParser(description="Summarize BAM read counts by gene biotype.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("gtf_file", help="Input GTF annotation file")
    parser.add_argument("-o", "--output", help="Output CSV file (defaults to BAM basename)")

    args = parser.parse_args()
    output_file = args.output or f"{os.path.splitext(args.bam_file)[0]}_biotype_counts.csv"

    # Process annotation
    genes = parse_gene(args.gtf_file)
    counts = count_biotype(genes, args.bam_file)
    df = pd.DataFrame(list(counts.items()), 
        columns=["Biotype", "Total Reads"])

    print(df)

if __name__ == "__main__":
    main()