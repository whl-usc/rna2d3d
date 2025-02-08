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
        temp_bam = f"{basename}.bam"
        subprocess.run(
            ["samtools", "view", "-@", str(n_threads), 
            "-bS", "-o", temp_bam, bam_path],
            capture_output=True, text=True, check=True
        )
        bam_path = temp_bam

    # Check if the sorted and indexed BAM file already exists
    sorted_bam = bam_path.replace(".bam", ".sorted.bam")
    index_file = sorted_bam + ".bai"

    if not os.path.exists(index_file):
        print(f"Index file not found. Sorting and indexing...")
        pysam.sort("-@", str(n_threads), "-o", sorted_bam, bam_path)
        pysam.index(sorted_bam)
        print(f"Input file sorted and indexed...")
    else:
        print(f"Sorted and indexed BAM file already exists: {sorted_bam}")

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
        print(f"Compressed file saved as: {file}.gz")
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
            4. CIGAR string + MD tag (concatenated with `column_20_value`)
            5. Repeat count of (position, CIGAR string)
            6. Strand ("+" or "-")
        - Prints the total number of reads processed.

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
    Converts BAM depth to BedGraph format efficiently using multiprocessing.

    This function calculates the depth of coverage from a BAM file using 
    `samtools depth`and converts output to BedGraph format in parallel
    using multiprocessing for improved efficiency.

    Parameters:
        bam_path (str): Path to the input BAM file.
        outfile_path (str): Path to the output BedGraph file (no extension).
        n_threads (int): Number of threads to use for processing.

    Returns:
        str: Path to the generated BedGraph file.

    Notes:
        - The output file is saved with a `.bedgraph` extension.
        - Lines with a depth value of 0 are skipped.
    """
    bedgraph_file = f"{outfile_path}.bedgraph"

    print(f"Processing {bam_path}...")
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
        description="""
Lossy compression of SAM/BAM file to zipped BED, BEDPE, or BEDGRAPH. 

'input_type' must be specified as one of the following:

    cont:       Non-gapped reads. Use the 'cont' option to process normal 
                RNA-seq data. 

    crssant:    Output SAM (or bam) file from crssant. Combined gap1 and trans 
                reads are clustered into unique groups and given a tag in
                specifying DG/NG/TG type. Note that input file usually contains
                'cliques' in the filename.

    gap1/gapm/homoe:  One or more gaps per read on the same chromosome.

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
        n_threads = str(nproc.stdout.strip())
    except subprocess.CalledProcessError:
        try:
            n_threads = os.cpu_count()
        except Exception as e:
            print(f"Error in retrieving CPU count: {e}. "
                f"Defaulting to 1 thread.")
            n_threads = 1

    basename = os.path.basename(args.infile).split('.bam')[0]
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