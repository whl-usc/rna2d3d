#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_04_30
Python:     python3.10
Script:     rna_inter_vis.py

Fix pritrans reads not being able to be clustered using CRSSANT.
"""

################################################################################
# Define version
__version__ = "1.0.0"

# Version notes
__update_notes__ = """
2.0.0
    -   Added instances for generating mini genome, processing pritrans files
        and preparing the file for CRSSANT DG assembly.

1.0.0
    -   Initial commit for function to parse files.
"""

################################################################################
# Import packages
import argparse
import math
import os
import pysam
import subprocess
import sys
from collections import defaultdict

################################################################################
# Define sub-functions for processing
def generate_mini_genome(input_fasta, region1, region2, output_fasta, genome_dir, 
                        n_threads, star_path="STAR"):
    """
    Extract two regions from input genome fasta file and create a mini genome using pysam.
    
    Args:
        input_fasta: Path to input fasta file
        region1: First region to extract (format: chr:start-end)
        region2: Second region to extract (format: chr:start-end)
        output_fasta: Path to output combined fasta file
        genome_dir: Directory to store STAR index
        threads: Number of threads for STAR
        star_path: Path to STAR executable
        
    Returns:
        Tuple of (new_chrom_name, region1_start, region1_end, region2_start, region2_end)
    """
    # Parse regions
    chr1, coords1 = region1.split(':')
    start1, end1 = map(int, coords1.split('-'))
    chr2, coords2 = region2.split(':')
    start2, end2 = map(int, coords2.split('-'))
    
    # Create new chromosome name
    new_chrom = f"{chr1}_{start1}-{end1}_{chr2}_{start2}-{end2}"
    
    # Use pysam to extract regions
    try:
        with pysam.FastaFile(input_fasta) as fa:
            seq1 = fa.fetch(reference=chr1, start=start1-1, end=end1)  # pysam uses 0-based
            seq2 = fa.fetch(reference=chr2, start=start2-1, end=end2)
            
            region1_len = len(seq1)
            region2_len = len(seq2)
            
            with open(output_fasta, 'w') as out_f:
                out_f.write(f">{new_chrom}\n")
                out_f.write(seq1 + "\n")
                out_f.write("N" * 20 + "\n")
                out_f.write(seq2 + "\n")
                
    except Exception as e:
        print(f"[ERROR] Failed to extract regions with pysam: {e}")
        sys.exit(1)
        
    # Calculate total genome length
    genome_length = region1_len + 20 + region2_len
    print(f"[INFO] Combined mini-genome length: {genome_length} bases.")
    
    # Calculate genomeSAindexNbases
    log2 = math.log2(genome_length)
    genomeSAindexNbases = int((log2 / 2) - 1)
    genomeSAindexNbases = min(genomeSAindexNbases, 14)
    
    print(f"[INFO] genomeSAindexNbases computed: {genomeSAindexNbases}")
    
    # Build STAR index
    print(f"[INFO] Building STAR index at '{genome_dir}'...")
    os.makedirs(genome_dir, exist_ok=True)
    
    try:
        subprocess.run([
            star_path,
            '--runThreadN', str(n_threads),
            '--runMode', 'genomeGenerate',
            '--genomeDir', genome_dir,
            '--genomeFastaFiles', output_fasta,
            '--genomeSAindexNbases', str(genomeSAindexNbases)
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to build STAR index: {e}")
        sys.exit(1)
        
    print("[INFO] STAR genome index built.")
    
    # Return new coordinates information
    return (new_chrom, 
            1,  							# Region 1 start in mini genome
            region1_len,  					# Region 1 end in mini genome
            region1_len + 21,  				# Region 2 start in mini genome (after 20 Ns)
            region1_len + 20 + region2_len) # Region 2 end in mini genome

def collapse_pritrans(input_pritrans, output_pritrans, new_chrom, r1_start, r1_end, r2_start, r2_end):
    """
    Collapse pritrans reads into single gapped alignments for the mini genome.
    
    Args:
        input_pritrans: Path to input pritrans file
        output_pritrans: Path to output collapsed pritrans file
        new_chrom: New chromosome name from mini genome
        r1_start: Start position of region 1 in mini genome
        r1_end: End position of region 1 in mini genome
        r2_start: Start position of region 2 in mini genome
        r2_end: End position of region 2 in mini genome
    """
    # Group alignments by read name
    read_alignments = defaultdict(list)
    
    with open(input_pritrans, 'r') as f:
        for line in f:
            if not line.startswith('@'):
                fields = line.strip().split('\t')
                read_name = fields[0]
                read_alignments[read_name].append(fields)
    
    # Process each read's alignments
    with open(output_pritrans, 'w') as out_f:
        with open(input_pritrans, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    out_f.write(line)
                else:
                    break
        
        for read_name, alignments in read_alignments.items():
            if len(alignments) == 1:
                # Single alignment - just write it with new coordinates
                fields = alignments[0]
                chrom = fields[1]
                start = int(fields[2])
                cigar = fields[5]
                
                # Convert to mini genome coordinates
                if chrom == new_chrom.split('_')[0]:  # Assuming first part is original chrom
                    new_start = r1_start if start < r1_end else r2_start
                else:
                    new_start = r2_start if start < r2_end else r1_start
                
                # Update fields
                fields[1] = new_chrom
                fields[2] = str(new_start)
                out_f.write('\t'.join(fields) + '\n')
            else:
                # Multiple alignments - create gapped alignment
                # This is simplified - may need more sophisticated logic
                primary_aln = alignments[0]
                fields = primary_aln.copy()
                
                # Update to new chromosome
                fields[1] = new_chrom
                
                # Set start to beginning of first region
                fields[2] = str(r1_start)
                
                # Create a gapped CIGAR (simplified)
                # This would need to be calculated based on actual alignments
                fields[5] = f"{r1_end-r1_start+1}M20N{r2_end-r2_start+1}M"
                
                out_f.write('\t'.join(fields) + '\n')

def prepare_crssant(input_sam, output_bam, remove_intermediate=True):
    """
    Convert SAM to sorted BAM and index for CRSSANT.
    
    Args:
        input_sam: Path to input SAM file
        output_bam: Path to output sorted BAM file (without .bam extension)
        remove_intermediate: Whether to remove intermediate files
    """
    # Convert SAM to BAM
    unsorted_bam = f"{output_bam}.unsorted.bam"
    with pysam.AlignmentFile(input_sam, 'r') as in_sam:
        with pysam.AlignmentFile(unsorted_bam, 'wb', header=in_sam.header) as out_bam:
            for read in in_sam:
                out_bam.write(read)
    
    # Sort BAM
    sorted_bam = f"{output_bam}.bam"
    pysam.sort('-o', sorted_bam, unsorted_bam)
    
    # Index BAM
    pysam.index(sorted_bam)
    
    # Clean up
    if remove_intermediate:
        os.remove(unsorted_bam)
        os.remove(input_sam)

def parse_arguments():
    """
    Set up command line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Process RNA-RNA interaction data for CRSSANT clustering.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')
    
    # Mini-genome generation command
    mini_parser = subparsers.add_parser('mini_genome', help='Generate mini genome')
    mini_parser.add_argument('-i', '--input_fasta', required=True, help='Input genome FASTA')
    mini_parser.add_argument('-r1', '--region1', required=True, help='First region (chr:start-end)')
    mini_parser.add_argument('-r2', '--region2', required=True, help='Second region (chr:start-end)')
    mini_parser.add_argument('-o', '--output_fasta', required=True, help='Output combined FASTA')
    mini_parser.add_argument('-g', '--genome_dir', required=True, help='STAR genome directory')
    mini_parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
    mini_parser.add_argument('--star_path', default='STAR', help='Path to STAR executable')
    
    # Pritrans collapsing command
    collapse_parser = subparsers.add_parser('collapse', help='Collapse pritrans reads')
    collapse_parser.add_argument('-i', '--input_pritrans', required=True, help='Input pritrans file')
    collapse_parser.add_argument('-o', '--output_pritrans', required=True, help='Output collapsed file')
    collapse_parser.add_argument('-c', '--new_chrom', required=True, help='New chromosome name')
    collapse_parser.add_argument('--r1_start', type=int, required=True, help='Region 1 start in mini genome')
    collapse_parser.add_argument('--r1_end', type=int, required=True, help='Region 1 end in mini genome')
    collapse_parser.add_argument('--r2_start', type=int, required=True, help='Region 2 start in mini genome')
    collapse_parser.add_argument('--r2_end', type=int, required=True, help='Region 2 end in mini genome')
    
    # CRSSANT preparation command
    crssant_parser = subparsers.add_parser('crssant', help='Prepare files for CRSSANT')
    crssant_parser.add_argument('-i', '--input_sam', required=True, help='Input SAM file')
    crssant_parser.add_argument('-o', '--output_bam', required=True, help='Output BAM base name')
    crssant_parser.add_argument('--keep_intermediate', action='store_true', 
                              help='Keep intermediate files')
    
    return parser.parse_args()

################################################################################
# Execute main

def main():
    args = parse_arguments()
    
   try:
        nproc = subprocess.run("nproc", shell=True, check=True, 
            text=True, capture_output=True)
        n_threads = int(nproc.stdout.strip())
    except subprocess.CalledProcessError:
        try:
            n_threads = os.cpu_count()
        except Exception as e:
            print(f"Error in retrieving CPU count: {e}. "
                f"Defaulting to 1 thread.")
            n_threads = 1

    if args.command == 'mini_genome':
        result = generate_mini_genome(
            args.input_fasta, args.region1, args.region2, 
            args.output_fasta, args.genome_dir, args.threads, args.star_path
        )
        print("Mini genome created successfully.")
        print(f"New chromosome: {result[0]}")
        print(f"Region 1 in mini genome: {result[1]}-{result[2]}")
        print(f"Region 2 in mini genome: {result[3]}-{result[4]}")
        
    elif args.command == 'collapse':
        collapse_pritrans(
            args.input_pritrans, args.output_pritrans, 
            args.new_chrom, args.r1_start, args.r1_end, args.r2_start, args.r2_end
        )
        print(f"Pritrans file collapsed and saved to {args.output_pritrans}")
        
    elif args.command == 'crssant':
        prepare_crssant(
            args.input_sam, args.output_bam, 
            not args.keep_intermediate
        )
        print(f"CRSSANT files prepared: {args.output_bam}.bam and index")
        
    else:
        print("Please specify a valid command. Use -h for help.")

if __name__ == "__main__":
    main()