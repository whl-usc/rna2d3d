#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_07_08
Python:     python3.10
Script:     rnaxrna.py

Fix pritrans reads not being able to be clustered using CRSSANT.
"""

################################################################################
# Define version
__version__ = "2.0.0"

# Version notes
__update_notes__ = """
2.0.0
    -   Added instances for generating mini genome, processing pritrans files
        and preparing the file for CRSSANT DG assembly.

1.0.0
    -   Set up new functions and logic.
    -   Include function to check for STAR aligner to be installed.
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

check_star()

def mini_genome(input_bam, ref_fasta, genome_dir, n_threads):
    """
    Extracts regions from input_bam and generates a mini-genome.
    
    Args:
        input_bam: Path to the input bam file to be processed
        ref_fasta: Path to the reference fasta to extract sequences from
        genome_dir: Directory to store STAR index, defaults to current directory
        n_threads: Number of threads for STAR 
    """
    os.makedirs(genome_dir, exist_ok=True)

    # Define start and end regions for the minigenome based on fastq
    if input_bam.endswith(".sam"):
        mode = "r"
    elif input_bam.endswith(".bam"):
        mode = "rb"
    else:
        raise ValueError("Unsupported file type: must be .sam or .bam")
    
    chroms, starts, ends = [], [], []

    with pysam.AlignmentFile(input_bam, mode) as bam:
        for read in bam:
            if read.is_unmapped:
                continue

            # Process records for the first RNA
            chroms.append(read.reference_name)
            starts.append(read.reference_start)
            ends.append(read.reference_end)

            # Handle records for second RNA...SA tagged.
            if read.has_tag("SA"):
                sa_fields = read.get_tag("SA").split(",")
                sa_pos = int(sa_fields[1]) - 1
                chroms.append(sa_fields[0])
                starts.append(int(sa_fields[1]))
                cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', sa_fields[3])
                end_len = sum(int(length) for length, op in cigar_ops if 
                    op in ('M', 'D', '=', 'X'))
                ends.append(sa_pos + end_len)

    regions = {}
    for chrom in set(chroms):
        chrom_starts = [s for c, s in zip(chroms, starts) if c == chrom]
        chrom_ends = [e for c, e in zip(chroms, ends) if c == chrom]
        regions[chrom] = (min(chrom_starts), max(chrom_ends))

    print(f"[INFO] Mini genome regions extracted from pritrans file...")
    print(f"{regions}")

    # Extract sequences from the reference fasta file
    fasta = pysam.FastaFile(ref_fasta)
    combined_seq = ""

    for chrom, (start, end) in sorted(regions.items()):
        seq = fasta.fetch(chrom, start, end)
        combined_seq += seq + "N" * 50
    combined_seq = combined_seq.rstrip("N")

    mini_genome_name = "mini_"+"_".join(k for k in regions.keys())
    with open(f"{genome_dir}/{mini_genome_name}.fa", "w") as out:
        out.write(f">{mini_genome_name}\n")
        out.write(f"{combined_seq}\n")

    # Calculate total genome length to calculate for genomeSAindexNbases
    genome_length = sum(end - start for start, end in regions.values())    
    print(f"[INFO] Combined mini-genome length: {genome_length} bases.")
    
    log2 = math.log2(genome_length)
    genomeSAindexNbases = int((log2 / 2) - 1)
    genomeSAindexNbases = min(genomeSAindexNbases, 14)   
    print(f"[INFO] genomeSAindexNbases: {genomeSAindexNbases}")

    # Build STAR index
    print(f"[INFO] Building STAR index at '{genome_dir}'...")    
    try:
        subprocess.run([
            'STAR',
            '--runThreadN', str(n_threads),
            '--runMode', 'genomeGenerate',
            '--genomeDir', genome_dir,
            '--genomeFastaFiles', str(f"{genome_dir}/{mini_genome_name}.fa"),
            '--genomeSAindexNbases', str(genomeSAindexNbases)
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to build STAR index: {e}")
        sys.exit(1)
        
    print("[INFO] STAR genome index built.")
    shutil.move("Log.out", genome_dir)

# mini_genome("SNORD133_tRNA_pri_crssant.sam", "hg38mask14add.fa", "./mini_genome", 6)
    
def bam_to_fastq(input_bam, output_fastq):
    with pysam.AlignmentFile(input_bam, "rb") as bam, open(output_fastq, "w") as fq:
        for read in bam:
            if read.is_unmapped or not read.is_secondary:
                fq.write(f"@{read.query_name}\n")
                fq.write(f"{read.query_sequence}\n")
                fq.write("+\n")
                fq.write(f"{read.qual}\n")

# bam_to_fastq("SNORD133_tRNA_pri_crssant.bam", "test.fastq")

def run_star(genome_dir, n_threads=6):
    fastq = "test.fastq"
    outprefix = "test"
    # Use new STAR
    print(f"[INFO] Aligning reads using STAR")    
    try:
        subprocess.run([
            'STAR',
            '--runThreadN', str(n_threads),
            '--runMode', 'alignReads',
            '--genomeDir', genome_dir,
            '--readFilesIn', fastq,
            '--outFileNamePrefix', outprefix,
            '--genomeLoad', 'NoSharedMemory',
            '--outReadsUnmapped', 'Fastx',
            '--outFilterMultimapNmax', '10',
            '--outFilterScoreMinOverLread', '0',
            '--outSAMattributes', 'All',
            '--outSAMtype', 'BAM', 'Unsorted', 'SortedByCoordinate',
            '--alignIntronMin', '1',
            '--scoreGap', '0',
            '--scoreGapNoncan', '0',
            '--scoreGapGCAG', '0',
            '--scoreGapATAC', '0',
            '--scoreGenomicLengthLog2scale', '-1',
            '--chimOutType', 'WithinBAM', 'HardClip',
            '--chimSegmentMin', '5',
            '--chimJunctionOverhangMin', '5',
            '--chimScoreJunctionNonGTAG', '0',
            '--chimScoreDropMax', '80',
            '--chimNonchimScoreDropMin', '20'
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to align using STAR: {e}")
        sys.exit(1)
        
    print("[INFO] STAR alignment completed.")

run_star("./mini_genome", "6")

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