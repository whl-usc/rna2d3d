"""
Contact:    wlee9829@gmail.com
Date:       2025_01_28
Python:     python3.10
Script:     bam_compressor.py

Converts a BAM file to modified BED/Bedgraph format while retaining
variation information. Merges records with identical start/stop locations.
"""

# Define version
__version__ = "1.1.0"

# Version notes
__update_notes__ = """
1.1.0
    -   Simplified processing function for gap1/gapm/homo file types to BED6 format.
1.0.0
    -   Initial commit, set up overall logic and functions.
"""

# Import packages
from collections import defaultdict
import argparse
import os
import pysam
import re

# Define functions
def check_and_index_bam(bam_path):
    """
    Ensure the BAM file is sorted and indexed. Create an index if it doesn't exist.
    """
    if not os.path.exists((bam_path).replace(".bam", "") + ".sorted.bam.bai"):
        print(f"Index file for {bam_path} not found. Sorting and indexing...")
        sorted_bam = bam_path.replace(".bam", ".sorted.bam")
        pysam.sort("-o", sorted_bam, bam_path)
        pysam.index(sorted_bam)

    return bam_path

def process_gap1_gapm_homo(bam_path, outfile_path):
    """
    Process gap1 files and extract relevant information.
    """
    repeat_counts = defaultdict(int)  # Count repeats of (position, cigar_items)
    with open(outfile_path+'.bed', 'w') as out_file:    
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue

                # SHORTER VERSION FOR BED6 OUTPUT
                chrom = bam.get_reference_name(read.reference_id)
                position = read.reference_start #chromStart in BED is 0-based, chromEnd is 1-based.
                cigar_items = read.cigarstring
                md_tag = read.get_tag("MD") if read.has_tag("MD") else "NA"
                strand = '-' if read.is_reverse else '+'

                # Increment repeat count for (position, cigar_items) and get count
                repeat_counts[(position, cigar_items)] += 1
                repeats = repeat_counts[(position, cigar_items)]

                out_file.write(f"{chrom}\t{position}\t.\t{read.cigarstring}_{md_tag}\t{repeats}\t{strand}\n")

def process_trans(bam_path, outfile_path):
    """
    Process gapm files and extract relevant information.
    """
    repeat_counts = defaultdict(int)
    with open(outfile_path+'.bed', 'w') as out_file:    
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue

def process_crssant(bam_path, outfile_path):
    repeat_counts = defaultdict(int)
    with open(outfile_path+'.bed', 'w') as out_file:    
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue

        # Extract the DG/NG/TG value
        column_20_value = None
        tags = read.tags  # Retrieve all optional tags
        if len(tags) >= 20:
            column_20_value = "_"+tags[19][1]
        elif len(tags) < 20:
            column_20_value = ""
            print(f"WARNING: NO DG/NG/TG TAG FOUND. CHECK COLUMN 20 FOR <--->")
            
        out_file.write(f"{chrom}\t{position}\t.\t{read.cigarstring}_{md_tag}{column_20_value}\t{repeats}\t{strand}\n")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Lossy compression of BAM data into modified BED format.")
    parser.add_argument("bam", help="Input BAM file.")
    parser.add_argument("type", help="gap1/gapm/trans/homo/cont/crssant")
    parser.add_argument("-C", "--clean", action="store_true", help="Removes the sorted BAM and index.")
    return parser.parse_args()

def main():
    args = parse_args()
    
    bam_path = check_and_index_bam(args.bam)
    basename = os.path.basename(args.bam).split('.')[0]
    type = str(args.type)
        
    if type in ("gap1", "gapm", "homo"):
        process_gap1_gapm_homo(bam_path, basename)
    elif type == "trans":
        pass
    elif type == "cont":
        pass
    elif type == "crssant"

    if args.clean:
        print("Cleaning up sorted BAM and index files...")
        os.remove(f"{basename}.sorted.bam")
        os.remove(f"{basename}.sorted.bam.bai")

if __name__ == "__main__":
    main()