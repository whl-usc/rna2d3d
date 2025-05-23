#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_01_28
Python:     python3.10
Script:     bam2bedz.py

Converts input BAM files to modified BED6/BEDPE/bedgraph format, where
count of "repeats" is a merged record with identical start/stop locations.

BED6: chromosome, position, ., cigarstring_md-tag, repeats, strand
BEDPE: chr1, start1, chr2, start2, cigar1_md1-cigar2_md2, 
        repeats, strand1, strand2 
BEDGRAPH: chromosome, start, end, count
"""

################################################################################
# Define version
__version__ = "2.3.0"

# Version notes
__update_notes__ = """
2.3.0
    -   Re-write input reading to accomodate for SAM as input.

2.2.0
    -   Added new function for continuous or "normal" read files.
    -   Overhaul of the BAM to bedgraph using multiprocessing 

2.1.1
    -   Optional flag to remove the sorted.bam and associated index file.

2.1.0
    -   Argument for analysis type made mandatory.

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
import os
import pysam
import re
import subprocess
import multiprocessing

################################################################################
# Define sub-functions for processing

def check_index_bam(bam_path, n_threads, basename):
    """
    Ensures the input file is in BAM format, sorted, and indexed. If input file 
    is SAM, convert to BAM before performing the sort and index.

    Parameters:
        bam_path (str): Path to the input file (BAM or SAM format).
        n_threads (int): Number of threads to use for sorting and indexing.
        basename (str): Base name from input file.

    Returns:
        str: Path to the sorted and indexed BAM file.

    Notes:
        - If a sorted and indexed BAM file already exists, no further
            processing is done.
    """
    if bam_path.endswith(".sam"):
        print(f"Provided file is in SAM format. Converting to BAM...")
        temp_bam = bam_path.replace(".sam", ".bam")
        subprocess.run(
            ["samtools", "view", "-@", str(n_threads), 
            "-bS", "-o", temp_bam, bam_path],
            capture_output=True, text=True, check=True
        )
        bam_path = temp_bam

    # Check if the sorted and indexed BAM file already exists
    sorted_bam = bam_path.replace(".bam", ".sorted.bam")
    index_file = sorted_bam + ".bai"

    if not os.path.exists(sorted_bam):
        print(f"Sorted BAM file not found. Sorting...")
        pysam.sort("-@", str(n_threads), "-o", sorted_bam, bam_path)
    if not os.path.exists(index_file):
        print(f"Index file not found. Indexing...")
        pysam.index(sorted_bam)
    else:
        print(f"Sorted and indexed BAM file already exists:\t{sorted_bam}\n")

    # Clean up the temporary BAM file if it was created
    if bam_path.endswith(".tmp.bam"):
        os.remove(bam_path)

    return sorted_bam


def compress_output(file, n_threads):
    """
    Compresses a file using bgzip with multi-threading support.

    Parameters:
        file (str): Path to the input file to be compressed.
        n_threads (int): Number of threads to use for compression.

    Notes:
        - The compressed file is saved with a `.gz` extension.
        - Requires `bgzip` from the htslib package to be installed.
    """
    try:
        subprocess.run(["bgzip", "-f", "-@", str(n_threads), file], check=True)
        print(f"Compressed file saved as:\t{file}.gz")
    except FileNotFoundError:
        print("Error: bgzip not found. Please install htslib.")
    except subprocess.CalledProcessError as e:
        print(f"Error: bgzip compression failed. {e}")


def compress_gap1_gapm_homo(bam_path, outfile_path, n_threads):
    """
    Processes a BAM file to extract relevant information and outputs a BED file.

    Returns:
        - A BED file ('outfile_path'.bed) containing:
            1. Chromosome
            2. Start position
            3. Placeholder (".")
            4. CIGAR string + MD tag
            5. Repeat count of (position, CIGAR string)
            6. Strand ("+" or "-")
        - Prints the total number of reads processed.

    Notes:
        - Unmapped reads are skipped.
        - Function tracks repeats for each (position, CIGAR string) pair.
        - Only the (position, CIGAR string) pairs with the highest repeat count are kept.
    """
    repeat_counts = defaultdict(int)  # Stores repeat counts per (position, CIGAR)
    read_data = {}  # Stores details for writing output
    total_reads = 0
    bed_file = f"{outfile_path}.bed"

    print(f"Processing {bam_path}...")

    # First pass: Count occurrences of (position, CIGAR)
    with pysam.AlignmentFile(bam_path, "rb", threads=n_threads) as bam:
        for read in bam:
            if read.is_unmapped:
                continue

            chrom = bam.get_reference_name(read.reference_id)
            position = read.reference_start
            cigar_items = read.cigarstring
            md_tag = read.get_tag("MD") if read.has_tag("MD") else "NA"
            strand = '-' if read.is_reverse else '+'

            key = (chrom, position, cigar_items)
            repeat_counts[key] += 1  # Increment repeat count
            total_reads += 1

            # Store BED entry if this is the highest count seen so far
            if repeat_counts[key] > read_data.get(key, (0, ""))[0]:
                read_data[key] = (
                    repeat_counts[key],  # Store the count
                    f"{chrom}\t{position}\t.\t{cigar_items}_{md_tag}\t{repeat_counts[key]}\t{strand}\n"
                )

    # Write only entries with max repeat counts
    with open(bed_file, 'w') as out_file:
        for _, (count, line) in read_data.items():
            out_file.write(line)

    print(f"Reads processed:\t{total_reads}")
    print(f"Results saved to:\t{bed_file}")

    return bed_file

def compress_trans(bam_path, outfile_path, n_threads):
    """
    Processes a BAM file to extract paired reads and outputs a BEDPE file.

    Returns:
        - A BEDPE file ('outfile_path'.bedpe) containing:
            1. Chromosome of read 1
            2. Start position of read 1
            3. Placeholder (".")
            4. Chromosome of read 2
            5. Start position of read 2
            6. Placeholder (".")
            7. CIGAR string + MD tag for read 1 and read 2 (joined with `_`)
            8. Repeat count of (chrom1, start1, cigar1, chrom2, start2, cigar2)
            9. Strand of read 1 ("+" or "-")
            10. Strand of read 2 ("+" or "-")
        - Prints the total number of paired reads processed.

    Notes:
        - Unmapped reads are skipped.
        - Paired reads are ordered by their start positions.
    """
    chimeric_reads = defaultdict(list)
    bedpe_lines = []
    chim_id = 0

    with pysam.AlignmentFile(bam_path, "rb", threads=n_threads) as bam:
        for read in bam:
            chimeric_reads[read.query_name].append(read)

    for read_name, reads in chimeric_reads.items():
        if len(reads) != 2:
            continue  # Only handle clean 2-part mappings

        r_primary = next(r for r in reads if not r.is_supplementary)
        r_supp = next(r for r in reads if r.is_supplementary)
        r1, r2 = r_primary, r_supp

        chrom1, start1 = r1.reference_name, r1.reference_start
        chrom2, start2 = r2.reference_name, r2.reference_start
        cigar1, cigar2 = r1.cigarstring, r2.cigarstring
        md1 = r1.get_tag("MD") if r1.has_tag("MD") else "NA"
        md2 = r2.get_tag("MD") if r2.has_tag("MD") else "NA"
        strand1 = '-' if r1.is_reverse else '+'
        strand2 = '-' if r2.is_reverse else '+'

        bedpe_line = (
            f"{chrom1}\t{start1}\t.\t"
            f"{chrom2}\t{start2}\t.\t"
            f"{cigar1}_{md1}-{cigar2}_{md2}\t"
            f"{chim_id}\t"
            f"{strand1}\t{strand2}\n"
        )
        bedpe_lines.append(bedpe_line)
        chim_id += 1

    out_path = f"{outfile_path}.bedpe"
    with open(out_path, 'w') as f:
        f.writelines(bedpe_lines)

    print(f"Chimeric (primary + supp) pairs written: {chim_id}")
    print(f"Output: {out_path}")
    return out_path


def process_line(line):
    """
    Convert a single samtools depth line to BedGraph format.
    """
    fields = line.strip().split()
    return (f"{fields[0]}\t{int(fields[1])-1}\t{fields[1]}\t{fields[2]}\n" 
        if len(fields) == 3 and int(fields[2]) > 0 else None)


def process_cont(bam_path, outfile_path, n_threads):
    """
    Converts BAM depth to BedGraph format efficiently using multiprocessing.

    This function calculates the depth of coverage from a BAM file using 
    samtools depth and converts output to BedGraph format in parallel
    using multiprocessing for improved efficiency.

    Parameters:
        bam_path (str): Path to the input BAM file.
        outfile_path (str): Path to the output BedGraph file (no extension).
        n_threads (int): Number of threads to use for processing.

    Returns:
        str: Path to the generated BedGraph file.

    Notes:
        - The output file is saved with a .bedgraph extension.
        - Lines with a depth value of 0 are skipped.
    """
    bedgraph_file = f"{outfile_path}.bedgraph"

    print(f"Processing {bam_path}...")
    result = subprocess.run(
        ["samtools", "view", "-@", str(n_threads), "-c", bam_path],
        capture_output=True, text=True, check=True
    )
    total_reads = int(result.stdout.strip())
    print(f"Reads processed:\t\t{total_reads}")
    
    with open(bedgraph_file, "w") as out_file, \
         subprocess.Popen(
            ["samtools", "depth", "-a", "-@", str(n_threads), bam_path], 
            stdout=subprocess.PIPE, text=True) as proc, \
         multiprocessing.Pool(int(n_threads)) as pool:
        
        for result in pool.imap_unordered(process_line, proc.stdout, 
            chunksize=500):
            if result:
                out_file.write(result)

    print(f"Results saved to:\t{bedgraph_file}")
    return bedgraph_file


def process_crssant(bam_path, outfile_path, n_threads):
    """π
    Processes a BAM file to extract relevant information and outputs a BED file.

    Returns:
        - A BED file ('outfile_path'.bed) containing:
            1. Chromosome
            2. Start position
            3. DG tag
            4. CIGAR string + MD tag
            5. Repeat count of (position, CIGAR string)
            6. Strand ("+" or "-")
        - Prints the total number of reads processed.

    Notes:
        - Unmapped reads are skipped.
        - Function tracks repeats for each (position, CIGAR string) pair.
    """
    repeat_counts = defaultdict(int)  # Track max repeat counts per (position, CIGAR)
    read_data = {} # Store corresponding BED entry for max repeat counts
    total_reads = 0
    bed_file = f"{outfile_path}.bed"

    print(f"Processing {bam_path}...")

    with pysam.AlignmentFile(bam_path, "rb", threads=n_threads) as bam:
        for read in bam:
            if read.is_unmapped:
                continue

            # Extract alignment details
            chrom = bam.get_reference_name(read.reference_id)
            position = read.reference_start
            cigar_items = read.cigarstring
            md_tag = read.get_tag("MD") if read.has_tag("MD") else "NA"
            strand = '-' if read.is_reverse else '+'

            # Process DG tag safely
            if read.has_tag("DG"):
                dg_raw = read.get_tag("DG")
                dg_tag = (dg_raw.replace(",", "_").rsplit(",\t", 1)[0] + ":" + 
                          dg_raw.rsplit(",", 1)[-1])
            else:
                dg_tag = "NA"

            # Create a unique key for tracking repeats
            key = (chrom, position, cigar_items, strand)

            # Update repeat count
            repeat_counts[key] += 1
            total_reads += 1

            # Store BED entry if this is the highest count seen so far
            if repeat_counts[key] > read_data.get(key, (0, ""))[0]:
                read_data[key] = (
                    repeat_counts[key],  # Store the count
                    f"{chrom}\t{position}\t{dg_tag}\t"
                    f"{cigar_items}_{md_tag}\t"
                    f"{repeat_counts[key]}\t{strand}\n"
                )

    # Write only entries with max repeat counts
    with open(bed_file, 'w') as out_file:
        for _, (count, line) in read_data.items():
            out_file.write(line)

    print(f"Reads processed:\t{total_reads}")
    print(f"Results saved to:\t{bed_file}")

    return bed_file


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
Lossy compression of SAM/BAM file to zipped BED, BEDPE, or BEDGRAPH. 

'input_type' must be specified as one of the following:

    cont:       Non-gapped reads. Use the 'cont' option to process normal 
                RNA-seq data. 

    crssant:    Output SAM (or bam) file from crssant. Combined gap1 and trans 
                reads are clustered into unique groups and given a tag in
                specifying DG/NG/TG type. Note that input file usually contains
                'cliques' in the filename.

    gap1/gapm/homo:  One or more gaps per read on the same chromosome.

    trans:      Trans-chromosome read. Individual arms are split into a paired  
                set of reads, spanning two lines. 
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("input_type", 
        help="Specify one of the following types of input: \
        cont crssant gap1 gapm homo trans")
    parser.add_argument("infile", 
        help="Input SAM or BAM file for compresion.")
    parser.add_argument("-r", "--remove", action="store_true", 
        help="Removes the sorted BAM and index.")

    return parser.parse_args()


def main():
    """
    Main function to process BAM file based on the specified analysis type.

    1. Parse input arguments.
    2. Ensure SAM/BAM file is sorted and indexed.
    3. Process the SAM/BAM file based on the name_type 
        (cont, crssant, gap1, gapm, homo, trans).
    4. Compress the output BED/BEDPE file using `bgzip`.
    5. Optionally clean up intermediate sorted BAM files.
    """
    args = parse_args()
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
    
    basename = os.path.splitext(os.path.basename(args.infile))[0]
    bam_file = check_index_bam(args.infile, n_threads, basename)
    analysis_type = str(args.input_type).lower()

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
        output = process_crssant(bam_file, basename, n_threads)
        compress_output(output, n_threads)

    else:
        print(f"Error: Unknown file type: {type}")
        return

    if args.remove:
        print("\nCleaning up sorted BAM and index files...")
        os.remove(f"{basename}.sorted.bam")
        os.remove(f"{basename}.sorted.bam.bai")

if __name__ == "__main__":
    main()