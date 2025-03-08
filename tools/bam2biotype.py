#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_03_04
Python:     python3.10
Script:     bam2biotypes.py
"""

################################################################################
# Define version
__version__ = "1.1.0"

# Version notes
__update_notes__ = """
1.2.0
    -   Dynamic thread allocation based on CPU count
    -   Automatic BAM indexing wherever .bai file is missing
    -   Improved read filtering and error handling

1.1.0
    -   Improved performance with multiprocessing for BAM read counting
    -   Added BAM index check and optimized the read filtering
    -   Enhanced GTF parsing efficiency
    -   Outputs results to CSV

1.0.0
    -   Initial commit, using pandas for GTF parsing.
"""
################################################################################
# Import packages
from collections import defaultdict
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
import subprocess

################################################################################
# Define sub-functions for processing
def get_n_threads():
    """
    Dynamically determine the number of threads to use.
    """
    try:
        return int(subprocess.check_output("nproc", shell=True).strip())
    except subprocess.CalledProcessError:
        try:
            return os.cpu_count() or 1
        except Exception as e:
            print(f"Error retrieving CPU count: {e}. Defaulting to 1 thread.")
            return 1


def ensure_bam_index(bam_file):
    """
    Check for BAM index (.bai) and create one if missing.
    """
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            if not bam.has_index():
                print(f"Indexing {bam_file}...")
                pysam.index(bam_file)
    except Exception as e:
        print(f"Error checking/indexing BAM: {e}")


def parse_gtf(gtf_file):
    """
    Parses a GTF file into a DataFrame, extracting gene_id and gene_biotype.
    """
    print(f"Processing {gtf_file} for gene_ids and gene_biotype...")
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
        "chrom", "source", "feature", "start", "end", "score", 
        "strand", "frame", "attributes"])

    df = pd.read_csv(gtf_file, sep="\t", comment="#", names=col_names, dtype=str)

    df = df[df["feature"] == "gene"]
    df["chrom"] = df["chrom"].replace(nc_to_chr)
    df["gene_id"] = df["attributes"].str.extract(r'gene_id "([^"]+)"')
    df["gene_biotype"] = df["attributes"].str.extract(r'gene_biotype "([^"]+)"')
    df = df.assign(gene_biotype=df["gene_biotype"].fillna("unknown"))

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
    
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    gtf_biotype = df.groupby("gene_id").agg({
        "chrom": "first",
        "start": "min",
        "end": "max",
        "strand": lambda x: x.iloc[0] if x.nunique() == 1 else "conflict",
        "gene_biotype": lambda x: x.iloc[0] if x.nunique() == 1 else "multiple"
    }).reset_index()

    # gene_counts = len(gtf_biotype["gene_id"])
    # print(gene_counts)
    # biotype_counts = gtf_biotype["gene_biotype"].value_counts()
    # print(biotype_counts)

    return gtf_biotype


def read_overlaps_gene(read_start, read_end, gene_start, gene_end, threshold=0.5):
    """
    Check if a read overlaps a gene by at least the given threshold.
    """
    overlap_length = max(0, min(read_end, gene_end) - max(read_start, gene_start))
    read_length = read_end - read_start
    
    return read_length > 0 and (overlap_length / read_length) >= threshold



def process_bam_by_gene(df, bam_file):
    """
    Process the BAM file and count uniquely mapped reads per gene.
    Then aggregate counts by biotype.

    Args:
        df (pd.DataFrame): DataFrame with gene info (columns: chrom, start, end, gene_id, biotype)
        bam_file (str): Path to BAM file

    Returns:
        gene_counts (dict): {gene_id: read_count}
        biotype_counts (dict): {biotype: read_count}
    """
    ensure_bam_index(bam_file)  

    gene_counts = defaultdict(int)
    biotype_counts = defaultdict(int)

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        valid_chroms = set(bam.references)

        for _, row in df.iterrows():
            chrom, start, end, gene_id, biotype = row["chrom"], int(row["start"]), int(row["end"]), row["gene_id"], row["gene_biotype"]

            if chrom not in valid_chroms:
                continue 

            seen_reads = set()  # Track reads per gene

            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                    continue
                
                read_start, read_end = read.reference_start, read.reference_end

                if read.query_name in seen_reads:  # Prevent duplicate counting per gene
                    continue

                if read_overlaps_gene(read_start, read_end, start, end):  
                    seen_reads.add(read.query_name)
                    gene_counts[gene_id] += 1
                    biotype_counts[biotype] += 1  

    return gene_counts, biotype_counts


def plot_stacked_biotype_counts(biotype_counts):
    """
    Plots a stacked bar chart with a single bar representing read counts split by biotype.

    Args:
        biotype_counts (dict): Dictionary with biotypes as keys and read counts as values.
    """
    # Convert to DataFrame
    df = pd.DataFrame(list(biotype_counts.items()), columns=["Biotype", "Read Count"])
    df = df.sort_values(by="Read Count", ascending=False)

    # Normalize to proportions for better cross-sample comparison
    df["Proportion"] = df["Read Count"] / df["Read Count"].sum()

    # Define colormap
    colors = plt.cm.viridis(np.linspace(0, 1, len(df)))

    # Plot stacked bar
    fig, ax = plt.subplots(figsize=(6, 6))

    bottom_values = np.zeros(1)  # Track cumulative proportions for stacking
    for i, (biotype, proportion, color) in enumerate(zip(df["Biotype"], df["Proportion"], colors)):
        ax.barh(["Reads"], proportion, left=bottom_values, color=color, label=biotype)
        bottom_values += proportion  # Stack each biotype on top

    # Labels & formatting
    ax.set_xlabel("Proportion of Reads")
    ax.set_title("Biotype Distribution as Stacked Proportions")

    # Add legend
    ax.legend(title="Biotype", loc="center left", bbox_to_anchor=(1, 0.5))

    # Save as PNG with 360 dpi
    plt.tight_layout()
    plt.savefig("stacked_biotype_read_count.png", dpi=360, bbox_inches="tight")


def main():
    """Main function for processing BAM file and summarizing biotype counts."""
    parser = argparse.ArgumentParser(description="Summarize BAM read counts by gene biotype.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("gtf_file", help="Input GTF annotation file")
    parser.add_argument("-o", "--output", help="Output CSV file (default: counts.csv)", default="counts.csv")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads (default: auto-detect)", default=get_n_threads())
    
    args = parser.parse_args()

    print(f"Using {args.threads} threads.")
    gtf_df = parse_gtf(args.gtf_file)
    gene_counts, biotype_counts = process_bam_by_gene(gtf_df, args.bam_file)

    plot_stacked_biotype_counts(biotype_counts)    
    
    # Save biotype counts
    biotype_output_df = pd.DataFrame(biotype_counts.items(), columns=["gene_biotype", "read_count"])
    biotype_output_df.to_csv(args.output, index=False)
    print(f"Biotype results saved to {args.output}")

    # Save per-gene read counts
    gene_output_file = args.output.replace(".csv", "_gene_counts.csv")
    gene_output_df = pd.DataFrame(gene_counts.items(), columns=["gene_id", "read_count"])
    gene_output_df.to_csv(gene_output_file, index=False)
    print(f"Gene count results saved to {gene_output_file}")

if __name__ == "__main__":
    main()