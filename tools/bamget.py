"""
contact:    wlee9829@gmail.com
date:       2024_02_09
python:     python3.10
script:     bamget.py

This script is used to prepare a BED file for CRSSANT analysis based on
the minimum coverage per nucleotide and annotation provided. Files should be
organized according to the mapping shell script from the rna2d3d repository.
"""

# Import Packages
from datetime import datetime
import argparse
import glob
import numpy as np
import os
import pandas as pd
import pysam
import subprocess
import sys
import textwrap

###########################################################################
def timenow():
    """
    Returns the current timestamp as a string.

    Returns:
        str: Current timestamp in format 'YYY-MM-DD HH:MM:SS'.
    """
    return str(datetime.now())[:-7]

def check_depth(bam_file, min_coverage):
    """
    Function that checks a specified provided bam_file to determine the 
    positions with minimum coverage as defined by the user. 

    Uses samtools depth function to determine reads at each nucleotide
    position, then collects positions >= min_coverage. Defining num_cores
    should speed up the depth consideration, but this is dependent on 
    the system which the script is run on. 

    Args:
        bam_file: PATH to BAM file.
        min_coverage: Integer defining minimum number of reads/nt before 
                      the region is considered covered.
    Returns:
        coverage_positions (list): [chromosome, positions]
    """
    coverage_positions = []
    num_cores = (os.cpu_count())/2; threads = f"-@ {num_cores}"
    with subprocess.Popen(['samtools', 'depth', threads, bam_file], stdout=subprocess.PIPE, text=True) as process:
        for line in process.stdout:
            chrom, pos, cov = line.split()
            if int(cov) >= int(min_coverage):
                coverage_positions.append((chrom, int(pos)))

    return coverage_positions

def split_bam_file(input_bam, min_coverage, skip_chromosome):
    """
    Splits a provided BAM file by their chromosome in column 3 and 
    determines coverage_positions using check_depth(). Collects the 
    values into a dictionary for further processing.  Future iterations 
    should rewrite the function to parallelize BAM file handling. 
    Files can also be written into memory rather than into intermediate
    *.tmp and *_sorted.bam if they are not necessary. If min_coverage
    is consistent, these files can also be deleted after processing. 
    Otherwise, keep the *_sorted.bam for faster parsing (reduces I/O).  

    Args:
        input_bam: PATH to BAM file.
        min_coverage: Integer defining minimum number of reads/nt before
                      a region is considered covered.

    Returns:
        high_coverage_positions (dictionary): list of chromosomes, positions.
    """
    
    try:
        read_bam = subprocess.check_output(['samtools', 'view', input_bam], text=True)
        header = subprocess.check_output(['samtools', 'view', '-H', input_bam], text=True)
        unique_values = {line.split('\t')[2] for line in read_bam.split('\n') if line}
        chromosomes = list(unique_values)
        high_coverage_positions = {}

        if skip_chromosome:
            for chrom in skip_chromosome:
                if chrom not in chromosomes:
                    print(f"Error: Chromosome '{chrom}' is not a valid chromosome from the BAM file. Please check for typos.")
                    sys.exit(1) # Exit the script with an error code.

        print(f"Checking positions with {min_coverage} or more reads.")
        # Use check_depth on each chromosome_sorted.bam file, using existing files where possible. 
        num_cores = (os.cpu_count())/2; threads = f"-@ {num_cores}"
        for chromosome in chromosomes:
            sorted_bam_file=f"{chromosome}_sorted.bam"
            if os.path.exists(sorted_bam_file):
                print(f"Sorted BAM file '{sorted_bam_file}' already exists for '{chromosome}'. Skipping split and sort.")
            else: 
                with open(f"{chromosome}.tmp", "w") as output_file:
                    output_file.write(header)
                    for line in read_bam.split('\n'):
                        if line.strip() and line.split('\t')[2] == chromosome:
                            output_file.write(line + '\n')

                # Sort the chromosome.tmp file using samtools sort
                subprocess.run(['samtools', 'sort', threads, f"{chromosome}.tmp", '-o', f"{chromosome}_sorted.bam"])
                os.remove(f"{chromosome}.tmp")
                print(f"Input BAM file has been succesfully split into {chromosome}_sorted.bam.")
            high_coverage_positions[chromosome] = check_depth(sorted_bam_file, min_coverage)

        return high_coverage_positions

    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to split BAM file and determine depth. Reason: {e}")
        return {}

def collapse_gene_regions(annotation_file):
    """
    Prepares a "collapsed" gene region file by parsing through a modified GTF 
    annotation file. See the help section for more details about preparing the
    file from a gencode GTF file. Future iterations may include handling of
    other annotation files.

    Args:
        annotation_file: input annotation file with gene information.

    Returns:
        gene_regions (dict): gene names as keys, chromosome, start, stop, and strand information as values.        
    """
    valid_chromosomes = set([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"])

    # Check if the output BED file already exists
    if os.path.exists("collapsed_genes.bed"):
        try:
            print(f"\nThe file 'collapsed_genes.bed' already exists. Reading gene regions from file.\n")
            gene_regions = {}
            with open("collapsed_genes.bed", 'r') as f:
                for line in f:
                    parts = line.split("\t")
                    chromosome = parts[0].strip("'\"")
                    start = int(parts[1])
                    stop = int(parts[2])
                    gene_name = parts[3]
                    strand = parts[4]
                    gene_regions[gene_name] = (chromosome, start, stop, strand)
        except:
            print(f"Error: Failed to parse the collapsed_genes.bed file. Check to see if the columns are as defined in the help section.")

    else:
        gene_regions = {}
        with open(annotation_file, 'r') as f:
            for line in f:
                parts = line.split("\t")
                chromosome = parts[0].strip("'\"")
                if chromosome not in valid_chromosomes:
                    chromosome = chromosome[3:]                
                start = int(parts[1])
                stop = int(parts[2])
                strand = parts[3]
                gene_name = parts[4]
                if gene_name not in gene_regions:
                    gene_regions[gene_name] = (chromosome, start, stop, strand)
                else:
                    gene_regions[gene_name] = (chromosome, min(start, gene_regions[gene_name][1]), max(stop, gene_regions[gene_name][2]), strand)

            # Write gene regions to BED file
            with open("collapsed_genes.bed", 'w') as file:
                for gene_name, (chromosome, start, stop, strand) in gene_regions.items():
                    file.write(f"{chromosome}\t{start}\t{stop}\t{gene_name}\t{strand}\n")

    return gene_regions

def overlap_gene_regions(high_coverage_positions, gene_regions):
    """
    Considers genomic regions that are overlapped by highly covered
    positions. For positions in a genomic coordinate region that have
    more than 5 reads, the gene is considered overlapped. Future 
    iterations may include a more complex calculation based on coverage
    percentage over gene length. However, this would need to take into
    consideration the kind of crosslinking or sequencing technology used.

    Args:
        high_coverage_positions (dictionary): list of chromosomes, positions that have >= min_coverage reads.
        gene_regions (dict): gene names as keys, chromosome, start, stop, and strand information as values.        

    Returns:
        overlapping_genes (set): Set of genes that have high_coverage_positions overlap. 
    """
    print(f"Checking gene regions for any overlap of highly covered positions.")
    overlapping_genes = set()
    stripped_positions = {chrom.strip("'\""): positions for chrom, positions in high_coverage_positions.items()}

    for gene_name, (gene_chrom, start, stop, strand) in gene_regions.items():
        gene_chrom = gene_chrom.strip("'\"")
        if gene_chrom in stripped_positions:
            covered_positions = sum(1 for pos in stripped_positions[gene_chrom] if start <= pos[1] <= stop)
            if covered_positions >= 5: # Change this cutoff if necessary.
                overlapping_genes.add(gene_name)

    return overlapping_genes

def write_bed_file(overlapping_genes, gene_regions, output):
    """
    Combines the functions to write output into a file for downstream processing.

    Args:
        overlapping_genes (set): Set of genes that have high_coverage_positions overlap.
        gene_regions (dict): Dictionary of all genomic regions from the annotation file.
        output (str): PATH to the output file. 
    
    Returns:
        None
    """
    mode = 'w'
    if os.path.exists(output):
        mode = 'a'

    existing_genes = set()
    if mode =='a':
        with open(output, 'r') as f:
            for line in f:
                existing_genes.add(line.split('\t')[3])

    with open(output, mode) as f:
        for gene_name in overlapping_genes:
            if gene_name in existing_genes:
                continue
            chrom, start, stop, strand = gene_regions[gene_name]
            line = f"{chrom}\t{start}\t{stop}\t{gene_name}\t1000\t{strand}\n"
            f.write(line)
            existing_genes.add(gene_name)

def main():
    parser = argparse.ArgumentParser(
        prog="bamget.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
    This script should be run in the parent directory where *_pri_crssant.bam is located
    after the mapping script as defined from the rna2d3d pipeline. It will process the
    BAM file that contains gapped_1 and trans reads, determining highly covered positions
    as defined by the user, then compare it to a user-provided annotation file to generate
    BED file of overlapping genes. The file can then be used directly for CRSSANT analysis.

    NOTE: Arguments should only be provided in the following order:

    1. input_bam:           PATH of the *pri_crssant.bam file that is generated after 
                            the mapping.sh script is used. Alternatively, any BAM
                            file can be supplied to check for read depth and genomic 
                            coorinate overlap.

    2. min_coverage:        Integer value that defines the minimum number of reads
                            at the nucleotide position before it is considered a
                            high coverage position. This value can be changed according
                            to the dataset. Start with a lower number and increase to
                            reduce the number of genes with only partial read overlap.

    3. annotation_file:     The annotation file is a modified GTF file which can be 
                            generated using gencode annotation files. Pre-made
                            annotation files must be renamed "collapsed_genes.bed". 
                            If working with PARIS/SHARC pipeline data, see the note below.

                            zcat gencode.v45.basic.annotation.gtf.gz | 
                            awk -F'\\t '$3 == "gene" && $9 ~ /gene_type "protein_coding"/ &&
                            $9 !~ /gene_name "ENSG0000"/ {split($9, a, "\\""); print $1 "\\t"
                            $4 "\\t" $5 "\\t" a[6] "\\t" $7}' > annotation.bed

                            NOTE: If missing specific genes that are defined as a separate 'chromosome' 
                            from the PARIS/SHARC pipelines, include them as separate lines in the annotation_file 
                            and check to see if present in output BED files:

                            Example genes with separate 'chromosomes':
                            ---------------------------------
                            RN7SK   1       331     RN7SK   +
                            RN7SL   1       288     RN7SL   +
                            RNU7    1       63      RNU7    +
                            RNY     1       112     RNY1    +
                            RNY     213     314     RNY2    +
                            RNY     415     510     RNY3    +
                            RNY     611     694     RNY4    +
                            U13     1       120     U13     +
                            U14AB   1       92      U14A    +
                            U14AB   193     283     U14B    +
                            U17     1       207     U17     +
                            U3      1       217     U3      +
                            U8      1       136     U8      +
                            hs12S   1       954     12S     +
                            hs16S   1       1559    16S     +
                            hs45S   3654    5523    18S     +
                            hs45S   6600    6757    5.8S    +
                            hs45S   7924    12994   28S     +
                            hs5S    1       121     hs5S    +
                            hssnRNA 1       164     U1      +
                            hssnRNA 265     451     U2      +
                            hssnRNA 552     696     U4      +
                            hssnRNA 797     902     U6      +
                            hssnRNA 1003    1118    U5      +
                            hssnRNA 1219    1353    U11     +
                            hssnRNA 1454    1603    U12     +
                            hssnRNA 1704    1833    U4atac  +
                            hssnRNA 1934    2058    U6atac  +
    
    4.  output:             Any name str() for the output file. It is recommended to include
                            min_coverage as part of the name to specify the cutoff values. Supplying
                            the same output filename will append to existing lines.  

    5.  -s, --skip-chromosome:

                            Optional argument. Provide a list of space separated chromosomes 
                            to omit when splitting the BAM file.

    ###########################################################################
    ###########################################################################
    """),
    usage="\npython3 %(prog)s [-h] [-rm] input_bam min_coverage annotation_file output [-s] chromosome ")
    
    parser.add_argument('-rm', '--remove', action='store_true', help='Optional parameter to remove the *_sorted.bam files.')
    parser.add_argument('-s', '--skip-chromosome', nargs='*', help='Optional parameter to skip specified chromosomes during processing.')
    parser.add_argument('input_bam')
    parser.add_argument('min_coverage')
    parser.add_argument('annotation_file')
    parser.add_argument('output')
    args = parser.parse_args()

    print(f"Job started at {timenow()}.\n")
    high_coverage_positions = split_bam_file(args.input_bam, args.min_coverage, args.skip_chromosome)
    gene_regions = collapse_gene_regions(args.annotation_file)
    overlapping_genes = overlap_gene_regions(high_coverage_positions, gene_regions)
    write_bed_file(overlapping_genes, gene_regions, args.output)

    remove_files = args.remove
    if remove_files:
        sorted_bam_files = glob.glob("*_sorted.bam")
        os.remove(collapsed_genes.bed)
        if sorted_bam_files:
            for file in sorted_bam_files:
                os.remove(file)
            print(f"Removed intermediate bam files.\n")
        else:
            print(f"No sorted BAM files found to remove.\n")
            

    print(f"Job completed at {timenow()}.\n")

if __name__ == "__main__":
    main()