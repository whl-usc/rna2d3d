#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_01_28
Python:     python3.10
Script:     bam2bedz.py

Converts input BAM files to modified BED6/BEDPE/bedgraph format, where
count of "repeats" is a merged record with identical start/stop locations.

BED6: chromosome, position, ., cigarstring_md-tag, repeats, strand
BEDPE: 
BEDGRAPH: chromosome, start, end, count
"""

################################################################################
# Define version

__version__ = "2.0.0"

# Version notes

__update_notes__ = """
2.0.0
    -   Order of arguments adjusted.  
    -   Added a separate function for trans reads.

1.1.0
    -   Simplified function for gap1/gapm/homo file types to BED6 format.
    -   Adjust count; BED chromStart is 0-based, chromEnd is 1-based.

1.0.0
    -   Initial commit, set up overall logic and functions.
"""

################################################################################
# Import packages

from collections import defaultdict
import argparse
import concurrent.futures
import os
import pysam
import re
import subprocess
import multiprocessing

################################################################################
# Define sub-functions for processing

def check_index_bam(bam_path, n_threads):
    """
    Ensure a BAM file is sorted and indexed.

    If an index (.bai) file does not exist for the BAM file, function:
    1. Sorts the BAM file and saves it as `*.sorted.bam`.
    2. Creates an index (`*.sorted.bam.bai`) for the sorted BAM file.

    Parameters:
        bam_path (str): Path to the input BAM file.

    Returns:
        str: Path to the sorted BAM file.
    """
    if not os.path.exists((bam_path).replace(".bam", "") + ".sorted.bam.bai"):
        print(f"Index file for {bam_path} not found.")
        print(f"Sorting and indexing...")
        sorted_bam = bam_path.replace(".bam", ".sorted.bam")
        pysam.sort("-@", n_threads, "-o", sorted_bam, bam_path)
        pysam.index(sorted_bam)
        print(f"{bam_path} sorted and indexed.")
    else:
        sorted_bam = bam_path.replace(".bam", ".sorted.bam")

    return sorted_bam

def compress_output(file, n_threads):
    try:
        subprocess.run(["bgzip", "-f", "-@", str(n_threads), file], check=True)
        print(f"Compressed file saved as: {file}.gz")
    except FileNotFoundError:
        print("Error: bgzip not found. Please install htslib.")
    except subprocess.CalledProcessError as e:
        print(f"Error: bgzip compression failed. {e}")


def compress_gap1_gapm_homo(bam_path, outfile_path, n_threads):
    """
    Processes a BAM file to extract relevant information and outputs a BED file.

    This function reads alignments from the BAM file, extracts key information 
    (chromosome, position, CIGAR string, MD tag, strand, and repeat count), and 
    writes it to a BED file in a modified BED6 format. 

    Parameters:
        bam_path (str): Path to the input BAM file.
        outfile_path (str): Base name for the output BED file

    Returns:
        str: Path to the generated BED file.

    Outputs:
        - A BED file ('outfile_path'.bed) containing:
            1. Chromosome
            2. Start position
            3. Placeholder (".")
            4. CIGAR string + MD tag (concatenated with `column_20_value`)
            5. Repeat count of (position, CIGAR string)
            6. Strand ("+" or "-")
        - Prints the total number of reads processed.

    Raises:
        FileNotFoundError: If the BAM file is missing or cannot be opened.
        KeyError: If required BAM tags (e.g., "MD") are unexpectedly absent.
        ValueError: If `bam_path` is not a valid BAM file.

    Notes:
        - Unmapped reads are skipped.
        - Function tracks repeats for each (position, CIGAR string) pair.
    """
    repeat_reads = defaultdict(int)  # Count repeats of (position, cigar_items)
    total_reads = 0 
    bed_file = (f"{outfile_path}.bed")

    print(f"Processing {bam_path}...")
    with open(bed_file, 'w') as out_file:    
        with pysam.AlignmentFile(bam_path, "rb", threads=n_threads) as bam:
            for read in bam:
                if read.is_unmapped:
                    continue

                # Retain only essential information for BED6 format
                chrom = bam.get_reference_name(read.reference_id)
                position = read.reference_start
                cigar_items = read.cigarstring
                md_tag = read.get_tag("MD") if read.has_tag("MD") else "NA"
                strand = '-' if read.is_reverse else '+'

                # Increment repeat count for (position, cigar_items) get count
                repeat_reads[(position, cigar_items)] += 1
                repeats = repeat_reads[(position, cigar_items)]

                out_file.write(
                    f"{chrom}\t{position}\t.\t"
                    f"{read.cigarstring}{'_'+md_tag}\t"
                    f"{repeats}\t{strand}\n"
                )

                total_reads += 1 

    print(f"Reads processed:\t\t{total_reads}")
    print(f"Results saved to:\t\t{bed_file}")

    return bed_file


def compress_trans(bam_path, outfile_path, n_threads):
    repeat_reads = defaultdict(int)
    total_pairs = 0  
    paired_reads = {}
    bedpe_file = (f"{outfile_path}.bedpe")

    print(f"Processing {bam_path}...")
    with open(bedpe_file, 'w') as out_file:
        with pysam.AlignmentFile(bam_path, "rb", threads=n_threads) as bam:
            for read in bam:
                if read.is_unmapped:
                    continue

                read_name = read.query_name

                if read_name not in paired_reads:
                    paired_reads[read_name] = []
                paired_reads[read_name].append(read)

            for read_name, reads in paired_reads.items():
                if len(reads) < 2:
                    continue

                read1, read2 = reads
                if read1.reference_start > read2.reference_start:
                    read1, read2 = read2, read1  # Ensure order start position

                # Extract key fields
                chrom1, start1 = read1.reference_name, read1.reference_start
                chrom2, start2 = read2.reference_name, read2.reference_start

                cigar1 = read1.cigarstring
                cigar2 = read2.cigarstring
                md1 = read1.get_tag("MD") if read1.has_tag("MD") else "NA"
                md2 = read2.get_tag("MD") if read2.has_tag("MD") else "NA"

                strand1 = '-' if read1.is_reverse else '+'
                strand2 = '-' if read2.is_reverse else '+'

                # Increment repeat count
                repeat_reads[(
                    chrom1, start1, cigar1, 
                    chrom2, start2, cigar2)] += 1
                repeats = repeat_reads[(
                    chrom1, start1, cigar1, 
                    chrom2, start2, cigar2)]

                # Write in BEDPE format
                out_file.write(
                    f"{chrom1}\t{start1}\t.\t"
                    f"{chrom2}\t{start2}\t.\t"
                    f"{cigar1}_{md1}-{cigar2}_{md2}\t"
                    f"{repeats}\t"
                    f"{strand1}\t{strand2}\n"
                )

                total_pairs += 1

    print(f"Paired reads processed:\t\t{total_pairs}")
    print(f"Results saved to:\t\t{bedpe_file}")

    return bedpe_file

def process_line(line):
    """
    Convert a single samtools depth line to BedGraph format.
    """
    fields = line.strip().split()
    return (f"{fields[0]}\t{int(fields[1])-1}\t{fields[1]}\t{fields[2]}\n" 
        if len(fields) == 3 and int(fields[2]) > 0 else None)

def process_cont(bam_path, outfile_path, n_threads):
    """
    Convert BAM depth to BedGraph format efficiently using multiprocessing.
    """
    bedgraph_file = f"{outfile_path}.bedgraph"

    print(f"Processing {bam_path}...")
    # Estimate total reads
    result = subprocess.run(
        ["samtools", "view", "-@", str(n_threads), "-c", bam_path],
        capture_output=True, text=True, check=True
    )
    total_reads = int(result.stdout.strip())
    print(f"Total reads: {total_reads}")
    
    with open(bedgraph_file, "w") as out_file, \
         subprocess.Popen(
            ["samtools", "depth", "-a", "-@", str(n_threads), bam_path], 
            stdout=subprocess.PIPE, text=True) as proc, \
         multiprocessing.Pool(int(n_threads)) as pool:
        
        # Process lines in parallel and write results
        for result in pool.imap_unordered(process_line, proc.stdout, 
            chunksize=500):
            if result:
                out_file.write(result)

    print(f"Results saved to: {bedgraph_file}")
    return bedgraph_file


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
            print(f"WARNING: DG/NG/TG TAG NOT FOUND. CHECK FOR <MD_tag>")
            
        out_file.write(
            f"{chrom}\t{position}\t.\t"
            f"{read.cigarstring}_{md_tag}{column_20_value}\t"
            f"{repeats}\t{strand}\n"
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Lossy compression for BAM into modified BED or BEDPE.")
    parser.add_argument("type", help="gap1/gapm/trans/homo/cont/crssant")
    parser.add_argument("infile", 
        help="Input BAM (compress) or BED (decompress) file.")
    parser.add_argument("-r", "--remove", action="store_true", 
        help="Removes the sorted BAM and index.")
    return parser.parse_args()


def main():
    """
    Main function to process BAM file based on the specified analysis type.

    1. Parse input arguments.
    2. Ensure BAM file is sorted and indexed.
    3. Process the BAM file based on the name_type (gap1, gapm, homo,
        trans, cont, crssant).
    4. Optionally clean up intermediate sorted BAM files.
    5. Optionally compress the output BED/BEDPE file using `bgzip`.
    """
    args = parse_args()
    try:
        # Attempt to get the number of processors using nproc
        nproc = subprocess.run("nproc", shell=True, check=True, 
            text=True, capture_output=True)
        n_threads = str(nproc.stdout.strip())
    except subprocess.CalledProcessError:
        try:
            n_threads = os.cpu_count()
        except Exception as e:
            print(f"Error in retrieving CPU count: {e}. "
                f"Defaulting to 1 thread.")
            n_threads = 1

    bam_file = check_index_bam(args.infile, n_threads)
    basename = os.path.basename(args.infile).split('.bam')[0]
    analysis_type = str(args.type)

    if analysis_type in ("gap1", "gapm", "homo"):
        output = compress_gap1_gapm_homo(bam_file, basename, n_threads)
        compress_output(output, n_threads)
    elif analysis_type == "trans":
        output = compress_trans(bam_file, basename, n_threads)
        compress_output(output, n_threads)
    elif analysis_type == "cont":
        output = process_cont(bam_file, basename, n_threads)
        compress_output(output, n_threads)
    elif analysis_type == "crssant":
        pass  # Placeholder for crssant analysis function
    else:
        print(f"Error: Unknown file type: {type}")
        return

    if args.remove:
        print("Cleaning up sorted BAM and index files...")
        os.remove(f"{basename}.sorted.bam")
        os.remove(f"{basename}.sorted.bam.bai")

if __name__ == "__main__":
    main()