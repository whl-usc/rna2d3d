#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_08_01
Python:     python3.10
Script:     genes2bed.py
"""
################################################################################
# Define version
__version__ = "2.0.0"

# Version notes
__update_notes__ = """
2.0.0
    -   Complete code refactor to use argparse and function blocks.
    -   Removed filter for protein_coding only genes.
"""
################################################################################

# Import packages
from datetime import datetime
import argparse
import os
import pandas as pd


def get_time():
    return str(datetime.now())[:-7]


def read_text(in_file, min_reads=1):
    """
    Read a tab-delimited file (e.g., ReadCount.txt) with a 'gene_Reads' column.
    Filter genes with reads >= min_reads and return the gene list.
    """
    print(f"{get_time()}\tReading file: {in_file} with min_reads={min_reads}")
    df = pd.read_csv(in_file, sep="\t")

    # Check required columns
    if "Gene" not in df.columns or "gene_Reads" not in df.columns:
        raise ValueError("Input file must contain 'Gene' and 'gene_Reads' columns")

    df_filtered = df[df["gene_Reads"] >= min_reads]
    gene_list = df_filtered["Gene"].tolist()

    print(f"{get_time()}\tFound {len(gene_list)} genes with >= {min_reads} reads")

    return gene_list


def read_anno(bed_file, genes_list, outprefix):
    """
    Reads a BED annotation file, filters rows for genes in genes_list,
    and writes a new BED file with collapsed regions per gene.
    """
    print(f"{get_time()}\tLoading annotation BED: {bed_file}")

    # Load BED file; if no header, specify column names explicitly
    bed = pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        names=["chromosome", "start", "stop", "strand", "gene", "type", "class"],
    )

    print(f"{get_time()}\tFiltering annotations for {len(genes_list)} genes")

    # Filter for genes in gene list
    bed_filtered = bed[bed["gene"].isin(genes_list)]

    # Collapse rows by gene to get min start and max stop per gene
    collapsed_rows = []
    for gene, group in bed_filtered.groupby("gene"):
        chromosome = group["chromosome"].iloc[0]
        strand = group["strand"].iloc[0]
        start = group["start"].min()
        stop = group["stop"].max()
        collapsed_rows.append(
            {
                "chromosome": chromosome,
                "start": start,
                "stop": stop,
                "gene": gene,
                "strand": strand,
            }
        )

    collapsed_df = pd.DataFrame(collapsed_rows)

    out_file = f"{outprefix}.bed"

    collapsed_df.to_csv(out_file, sep="\t", header=False, index=False)

    print(
        f"{get_time()}\tWrote collapsed BED with {len(collapsed_df)} genes to {out_file}"
    )

    return out_file


def parse_args():
    parser = argparse.ArgumentParser(
        description="""This script generates a collapsed gene BED file from
        a ReadCount.txt file output by CountCdsUtr.py and an annotation BED file. It filters genes by minimum reads."""
    )

    parser.add_argument("readcount", help="ReadCount.txt output from CountCdsUtr.py")
    parser.add_argument("annotation", help="Annotation BED file with gene regions")
    parser.add_argument(
        "min_reads",
        type=int,
        help="Minimum reads for gene inclusion (default=1)",
    )
    parser.add_argument("outprefix", help="Output file prefix")

    return parser.parse_args()


def main():
    args = parse_args()

    if "ReadCount.txt" not in args.readcount:
        print(
            "Warning: The readcount file does not seem to be from CountCdsUtr.py output"
        )

    genes = read_text(args.readcount, min_reads=args.min_reads)
    collapsed_bed = read_anno(args.annotation, genes, args.outprefix)


if __name__ == "__main__":
    main()
