#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_08_21
Python:     python3.10
Script:     rnaxrna.py

Script re-processes pritrans gapped files to be mapped to a single mini-genome.
1) Check for STAR and samtools/pysam, input files as FASTA/BAM
2) Process input files for miniGenome generation
3) Re-map BAM to mini-genome
4) Show instructions for running crssant_birch...print med seglen
"""

################################################################################
# Define version
__version__ = "2.0.0"

# Version notes
__update_notes__ = """
1.0.0
    -   Set up logic and workflow.
    -   Set up a function to check for STAR aligner installation.
"""

################################################################################
# Import packages
from collections import defaultdict
import argparse
import math
import os
import pysam
import re
import shutil
import subprocess
import sys


################################################################################
# Define sub-functions for processing
def check_star():
    """
    Checks if STAR aligner runs properly by calling --help.

    Returns:
        None if STAR works.
        Prints warning if STAR fails.
    """
    if shutil.which("STAR") is None:
        raise EnvironmentError("STAR not found in PATH")


# check_star()


def read_fasta(file, out_file=None):
    """
    Read a FASTA containing exactly two sequences and merge them into a
    single  entry, separated by 50 Ns.

    Inputs:
        FASTA with two sequences (e.g. from NCBI Nucleotide:
        https://www.ncbi.nlm.nih.gov/nuccore)

    Example:
        >Gene1
        AAAAAA
        >Gene2
        CCCCCC

    Becomes:
        >Gene1_Gene2
        AAAAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCCCCCC
    """
    genes = {}
    current_gene = None
    seq_chunks = []
    # Parse fasta into dict
    with open(file, "rt") as fasta:
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_gene:
                    genes[current_gene] = "".join(seq_chunks)

                header = line[1:].strip()
                match = re.search(r"\(([^)]+)\)", header)
                if match:
                    current_gene = match.group(1)
                else:
                    current_gene = header.split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip().upper())
        if current_gene:
            genes[current_gene] = "".join(seq_chunks)

    # Combine sequences with 50 Ns
    gene_names = list(genes.keys())
    if len(gene_names) != 2:
        raise ValueError("Need exactly two sequences in FASTA.")

    combined_name = "_".join(gene_names)
    combined_seq = ("N" * 50).join(genes[name] for name in gene_names)

    # Write output fasta
    if out_file:
        with open(out_file, "w") as out:
            out.write(f">{combined_name}\n")
            for i in range(0, len(combined_seq), 60):
                out.write(combined_seq[i : i + 60] + "\n")
        return out_file
    else:
        fasta_str = f">{combined_name}\n"
        for i in range(0, len(combined_seq), 60):
            fasta_str += combined_seq[i : i + 60] + "\n"
        # print(fasta_str)

        return fasta_str


read_fasta("./gapdh_18s.fasta")


def mini_genome(ref_fasta, genome_dir, n_threads=1):
    """
    This function takes a reference FASTA file and builds a STAR genome index.

    Args:
        ref_fasta (str): Reference FASTA used for genome generation function.
        genome_dir (str): Path where mini-genome will be created.
        n_threads (int): Number of threads used for STAR

    """


mini_genome("SNORD133_tRNA_pri_crssant.sam", "hg38mask14add.fa", "./mini_genome", 6)


# def run_star(genome_dir, n_threads=6):
#     fastq = "test.fastq"
#     outprefix = "test"
#     # Use new STAR
#     print(f"[INFO] Aligning reads using STAR")
#     try:
#         subprocess.run(
#             [
#                 "STAR",
#                 "--runThreadN",
#                 str(n_threads),
#                 "--runMode",
#                 "alignReads",
#                 "--genomeDir",
#                 genome_dir,
#                 "--readFilesIn",
#                 fastq,
#                 "--outFileNamePrefix",
#                 outprefix,
#                 "--genomeLoad",
#                 "NoSharedMemory",
#                 "--outReadsUnmapped",
#                 "Fastx",
#                 "--outFilterMultimapNmax",
#                 "10",
#                 "--outFilterScoreMinOverLread",
#                 "0",
#                 "--outSAMattributes",
#                 "All",
#                 "--outSAMtype",
#                 "BAM",
#                 "Unsorted",
#                 "SortedByCoordinate",
#                 "--alignIntronMin",
#                 "1",
#                 "--scoreGap",
#                 "0",
#                 "--scoreGapNoncan",
#                 "0",
#                 "--scoreGapGCAG",
#                 "0",
#                 "--scoreGapATAC",
#                 "0",
#                 "--scoreGenomicLengthLog2scale",
#                 "-1",
#                 "--chimOutType",
#                 "WithinBAM",
#                 "HardClip",
#                 "--chimSegmentMin",
#                 "5",
#                 "--chimJunctionOverhangMin",
#                 "5",
#                 "--chimScoreJunctionNonGTAG",
#                 "0",
#                 "--chimScoreDropMax",
#                 "80",
#                 "--chimNonchimScoreDropMin",
#                 "20",
#             ],
#             check=True,
#         )
#     except subprocess.CalledProcessError as e:
#         print(f"[ERROR] Failed to align using STAR: {e}")
#         sys.exit(1)

#     print("[INFO] STAR alignment completed.")


# run_star("./mini_genome", "6")

# def prepare_crssant(input_sam, remove_intermediate=True):
#     """
#     Convert SAM to sorted BAM and index for CRSSANT using pysam.

#     Args:
#         input_sam: Path to input SAM file
#         remove_intermediate: Whether to remove intermediate files (e.g. unsorted BAM, SAM)
#     """
#     # Derive output paths from input prefix
#     input_prefix = os.path.splitext(input_sam)[0]
#     unsorted_bam = f"{input_prefix}_unsorted.bam"
#     sorted_bam = f"{input_prefix}_sorted.bam"

#     # Convert SAM to BAM (unsorted)
#     with pysam.AlignmentFile(input_sam, "r") as in_sam:
#         with pysam.AlignmentFile(unsorted_bam, "wb", header=in_sam.header) as out_bam:
#             for read in in_sam:
#                 out_bam.write(read)

#     # Sort BAM
#     pysam.sort("-o", sorted_bam, unsorted_bam)

#     # Index BAM
#     pysam.index(sorted_bam)

#     # Remove intermediate files
#     if remove_intermediate:
#         os.remove(unsorted_bam)
#         os.remove(input_sam)

# def parse_arguments():
#     """
#     Set up command line arguments.
#     """
#     parser = argparse.ArgumentParser(
#         description='Process RNA-RNA interaction data for CRSSANT clustering.',
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#     subparsers = parser.add_subparsers(dest='command', help='Sub-commands')

#     # Mini-genome generation command
#     mini_parser = subparsers.add_parser('mini_genome', help='Generate mini genome')
#     mini_parser.add_argument('-i', '--input', required=True, help='Input genome FASTA')
#     mini_parser.add_argument('-r1', '--region1', required=True, help='First region (chr:start-end)')
#     mini_parser.add_argument('-r2', '--region2', required=True, help='Second region (chr:start-end)')
#     mini_parser.add_argument('-o', '--output', required=True, help='Output combined FASTA')
#     mini_parser.add_argument('-g', '--genome_dir', required=True, help='STAR genome directory')
#     mini_parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
#     mini_parser.add_argument('--star_path', default='STAR', help='Path to STAR executable')

#     # Pritrans collapsing command
#     collapse_parser = subparsers.add_parser('collapse', help='Collapse pritrans reads')
#     collapse_parser.add_argument('-i', '--input', required=True, help='Input pritrans file')
#     collapse_parser.add_argument('-o', '--output', required=True, help='Output collapsed file')
#     collapse_parser.add_argument('-c', '--new_chrom', required=True, help='New chromosome name')
#     collapse_parser.add_argument('--r1_start', type=int, required=True, help='Region 1 start in mini genome')
#     collapse_parser.add_argument('--r1_end', type=int, required=True, help='Region 1 end in mini genome')
#     collapse_parser.add_argument('--r2_start', type=int, required=True, help='Region 2 start in mini genome')
#     collapse_parser.add_argument('--r2_end', type=int, required=True, help='Region 2 end in mini genome')

#     # CRSSANT preparation command
#     crssant_parser = subparsers.add_parser('crssant', help='Prepare files for CRSSANT')
#     crssant_parser.add_argument('-i', '--input', required=True, help='Input SAM file')
#     crssant_parser.add_argument('-o', '--output', required=True, help='Output BAM base name')
#     crssant_parser.add_argument('--keep_intermediate', action='store_true',
#                               help='Keep intermediate files unsorted BAM files')

#     return parser.parse_args()

# ################################################################################
# # Execute main

# def main():
#     args = parse_arguments()

#     try:
#         nproc = subprocess.run("nproc", shell=True, check=True,
#             text=True, capture_output=True)
#         n_threads = int(nproc.stdout.strip())
#         print(f"Using {n_threads} CPUs for processing.")
#         if n_threads < 8:
#             print(f"WARNING: STAR requires >=16G RAM for alignment.")
#     except subprocess.CalledProcessError:
#         try:
#             n_threads = os.cpu_count()
#         except Exception as e:
#             print(f"Error in retrieving CPU count: {e}. Defaulting to 1 thread."
#                 f"WARNING: STAR takes longer to run with a single thread!!!")
#             n_threads = 1

#     check_star()

#     if args.command == 'mini_genome':
#         result = generate_mini_genome(
#             args.input_fasta, args.region1, args.region2,
#             args.output_fasta, args.genome_dir, args.threads, args.star_path
#         )
#         print("Mini genome created successfully.")
#         print(f"New chromosome: {result[0]}")
#         print(f"Region 1 in mini genome: {result[1]}-{result[2]}")
#         print(f"Region 2 in mini genome: {result[3]}-{result[4]}")

#     elif args.command == 'collapse':
#         collapse_pritrans(
#             args.input_pritrans, args.output_pritrans,
#             args.new_chrom, args.r1_start, args.r1_end, args.r2_start, args.r2_end
#         )
#         print(f"Pritrans file collapsed and saved to {args.output_pritrans}")

#     elif args.command == 'crssant':
#         prepare_crssant(
#             input_sam=args.input_sam,
#             output_bam=args.output_bam,
#             remove_intermediate=not args.keep_intermediate
#         )
#         print(f"CRSSANT files prepared: {args.output_bam}.bam and index")

#     else:
#         print("Please specify a valid command. Use -h for help.")

# if __name__ == "__main__":
#     main()
