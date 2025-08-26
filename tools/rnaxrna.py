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
__version__ = "1.0.0"

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
from datetime import datetime
import math
import numpy as np
import os
import pysam
import re
import shutil
import subprocess
import sys


################################################################################
# Define sub-functions for processing
def timenow():
    return f"{datetime.now().strftime('%b %d %H:%M:%S')} ---"


def filter_bam(input_bam, chrA, chrB, outdir=None):
    """
    Filter a BAM file to keep reads mapping between chrA and chrB.

    Args:
        input_bam (str): Input BAM file.
        chrA (str): First chromosome/contig.
        chrB (str): Second chromosome/contig.
        outdir (str, optional): Output directory. Defaults to same as input.

    Returns:
        str: Path to filtered BAM file.
    """
    outname = os.path.join(outdir, f"{chrA}_{chrB}_filtered.bam")

    with pysam.AlignmentFile(input_bam, "rb") as in_bam, pysam.AlignmentFile(
        outname, "wb", template=in_bam
    ) as out_bam:

        for read in in_bam:
            if read.is_unmapped:
                continue

            ref = in_bam.get_reference_name(read.reference_id)
            # For reads with SA tag (chimeric), it may contain other chromosomes
            sa = read.get_tag("SA") if read.has_tag("SA") else ""

            # Check if read maps between chrA and chrB in any orientation
            if (ref == chrA and chrB in sa) or (ref == chrB and chrA in sa):
                out_bam.write(read)

        print(f"{timenow()} Filtered BAM saved as: {os.path.basename(outname)}")

    return outname


def read_fasta(file, out_fasta, out_bed):
    """
    Read a FASTA with multiple sequences and merge into one with Ns between.
    Optionally create a BED file with coordinates of each original sequence.

    Args:
        file: Input FASTA
        out_fasta: Path to write combined FASTA
        out_bed: Path to write BED (gene_name, start, end)
        gap_len: Number of Ns between sequences

    Inputs:
        One FASTA file containing at least one entry from NCBI Nucleotide:
            https://www.ncbi.nlm.nih.gov/nuccore

        >NM_001256799.3 Homo sapiens glyceraldehyde-3-phosphate dehydrogenase (GAPDH), transcript variant 2, mRNA
        GTCCGGATGCTGCGCCTGCGGTAGAGCGGCCGCCATGTTGCAACCGGGAAGGAAATGAATGGGCAGCCGT

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

    gene_names = list(genes.keys())
    combined_name = "_".join(gene_names)
    combined_seq = ("N" * 20).join(genes[name] for name in gene_names)

    print(f"{timenow()} Preparing fasta genes: {combined_name}")
    with open(out_fasta, "w") as fasta:
        fasta.write(f">{combined_name}\n")
        for i in range(0, len(combined_seq), 60):
            fasta.write(combined_seq[i : i + 60] + "\n")

    with open(out_bed, "w") as bed:
        start = 0
        for name in gene_names:
            seq_len = len(genes[name])
            bed.write(f"{combined_name}\t{start}\t{start + seq_len}\t{name}\t1000\t+\n")
            start += seq_len + 50

    return out_fasta


def check_star():
    """
    Checks if STAR aligner is available and runnable.
    """
    STAR_path = shutil.which("STAR")
    if STAR_path is None:
        print(f"{timenow()} STAR not found in PATH; check installation.")
        sys.exit(1)

    try:
        subprocess.run(
            [STAR_path, "--help"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        print(f"{timenow()} Using STAR at {STAR_path}")
    except subprocess.CalledProcessError:
        print(f"{timenow()} Tried using STAR at {STAR_path}, but failed to run.")
        sys.exit(2)


def mini_genome(ref_fasta, genome_dir, n_threads=1):
    """
    Builds a STAR genome index from reference FASTA file.

    Args:
        ref_fasta (str): Reference FASTA used for genome generation function.
        genome_dir (str): Path where mini-genome will be created.
        n_threads (int): Number of threads used for STAR

    """
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)

    # Estimate --genomeSAindexNbases for small genomic indices
    genome_size = 0
    with open(ref_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            genome_size += len(line.strip())

    SAindexNbases = min(14, int(math.log2(genome_size) / 2 - 1))

    # fmt: off
    cmd = [
        "STAR", 
        "--runMode", "genomeGenerate", 
        "--genomeDir", genome_dir, 
        "--genomeFastaFiles", ref_fasta, 
        "--runThreadN", str(n_threads), 
        "--genomeSAindexNbases", str(SAindexNbases),
    ]  # fmt: on

    print(f"\nRunning STAR genomeGenerate.")
    subprocess.run(cmd, check=True)
    print(f"\n{timenow()} STAR genome index built at {genome_dir}\n")


def bam_to_gapped_fastq(bam_file, fastq_out, gap_len=10):
    """
    Collapse BAM reads into two-arm gapped FASTQ.
    If multiple supplementary alignments exist, only the two strongest arms are kept.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    out = open(fastq_out, "w")

    read_groups = {}

    # Group all alignments by qname
    for read in bam:
        if read.is_unmapped:
            continue
        qname = read.query_name
        read_groups.setdefault(qname, []).append(read)

    for qname, group in read_groups.items():
        if len(group) < 2:
            continue

        group.sort(key=lambda r: r.query_alignment_length, reverse=True)
        r1, r2 = group[:2]

        seq = r1.query_sequence + r2.query_sequence
        qual = "".join(chr(q + 33) for q in r1.query_qualities) + "".join(
            chr(q + 33) for q in r2.query_qualities
        )

        out.write(f"@{qname}\n{seq}\n+\n{qual}\n")

    bam.close()
    out.close()


def map_to_minigenome(combined_fastq, genome_dir, n_threads=8, out_prefix=None):
    """
    Map combined fastq to mini-genome with custom STAR parameters, and extract primary alignments.

    Args:
        input_bam (str): Input BAM file with paired reads.
        genome_dir (str): STAR genome index directory.
        n_threads (int): Threads for STAR.
        out_prefix (str): Prefix for STAR outputs.
    """

    # STAR mapping with custom parameters
    # fmt: off
    star_cmd = [
        "STAR",
        "--runThreadN", str(n_threads),
        "--runMode", "alignReads",
        "--genomeDir", genome_dir,
        "--readFilesIn", combined_fastq,
        "--outFileNamePrefix", f"{out_prefix}_",
        "--genomeLoad", "NoSharedMemory",
        "--outReadsUnmapped", "Fastx",
        "--outFilterMultimapNmax", "10",
        "--outFilterScoreMinOverLread", "0",
        "--outFilterMismatchNoverLmax", "1.0",
        "--outSAMattributes", "All",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--alignIntronMin", "1",
        "--alignSJDBoverhangMin", "1",
        "--alignSJoverhangMin", "1",
        "--scoreGap", "0",
        "--scoreGapNoncan", "0",
        "--scoreGapGCAG", "0",
        "--scoreGapATAC", "0",
        "--scoreGenomicLengthLog2scale", "-1",
        "--chimOutType", "WithinBAM", "HardClip",
        "--chimSegmentMin", "5",
        "--chimJunctionOverhangMin", "1",
        "--chimScoreJunctionNonGTAG", "0",
        "--chimScoreDropMax", "80",
        "--chimNonchimScoreDropMin", "20",
        "--chimScoreSeparation", "0"
    ]  # fmt: on

    print("Running STAR mapping to mini_genome.")
    subprocess.run(star_cmd, check=True)

    # Extract primary alignments
    star_bam = f"{out_prefix}_Aligned.sortedByCoord.out.bam"
    pri_bam = f"{out_prefix}_pri.bam"
    if os.path.exists(star_bam):
        with pysam.AlignmentFile(star_bam, "rb") as bam:
            with pysam.AlignmentFile(pri_bam, "wb", template=bam) as out_bam:
                for read in bam:
                    # Keep all header lines and reads with <21 fields (non-chimeric)
                    if read.is_secondary or read.is_supplementary:
                        continue
                    out_bam.write(read)
        print(f"\n{timenow()} Primary alignments saved as {pri_bam}")
    else:
        print(f"\n{timenow()} STAR BAM {star_bam} not found. Check STAR run.")


def cleanup(outdir, out_prefix):
    os.remove(f"{outdir}/{out_prefix}.fastq")
    os.remove(f"Log.out")
    for suffix in ["Log.out", "Log.progress.out", "SJ.out.tab", "Unmapped.out.mate1"]:
        os.remove(f"{outdir}/{out_prefix}_{suffix}")


def get_segment_lengths(bam_file):
    """
    Extract lengths of all mapped segments from a BAM file (or SAM).

    Args:
        bam_file (str): Input BAM or SAM file path.

    Returns:
        seg_lengths (list[int]): List of mapped segment lengths.
        median_len (float, optional): Median segment length (if return_median=True).
    """
    seg_lengths = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            # CIGAR is a list of tuples: (operation, length)
            # Operation codes: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 7/=, 8=X
            # Split by 'N' operations to get segments
            start = 0
            seg_len = 0
            for op, length in read.cigartuples:
                if op == 0 or op == 7 or op == 8:  # M/=X
                    seg_len += length
                elif op == 3:  # N
                    if seg_len > 0:
                        seg_lengths.append(seg_len)
                    seg_len = 0
            if seg_len > 0:
                seg_lengths.append(seg_len)

        seglen_median = float(np.median(seg_lengths)) if seg_lengths else 0
        # print(f"Alignments with seglen: {seg_lengths}")

    return seglen_median


def crssant_birch(bam_file, bed_file, outdir, n_threads, seglen_median):
    """
    Sets up and runs CRSSANT_birch.
    """

    # Derive output paths from input prefix
    print(f"{timenow()} Median seglen used for CRSSANT_birch: {seglen_median}")


def parse_arguments():
    """
    Set up command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Process RNA-RNA interaction data for CRSSANT clustering.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--fasta",
        "-f",
        help="Input FASTA file containing sequences to merge into a mini-genome.",
    )
    parser.add_argument(
        "--bam",
        "-b",
        required=True,
        help="Input BAM file containing trans-paired reads to convert into gapped FASTQ.",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        default="./rnaxrna_out",
        help="Output directory for mini-genome and mapping results.",
    )
    parser.add_argument(
        "--prefix", "-p", default="rnaxrna", help="Prefix for STAR output files."
    )
    parser.add_argument(
        "--gap-len",
        "-g",
        type=int,
        default=20,
        help="Number of Ns to insert between paired reads in recombined FASTQ.",
    )
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=None,
        help="Number of threads for STAR (auto-detected if not provided).",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep intermediate files (default: remove).",
    )
    parser.add_argument(
        "--filter",
        nargs=2,
        metavar=("CHR_A", "CHR_B"),
        help="Filter BAM to reads connecting CHR_A and CHR_B.",
    )
    parser.add_argument(
        "--version", "-v", action="version", version=f"%(prog)s {__version__}"
    )

    return parser.parse_args()


################################################################################
# Execute main
def main():
    args = parse_arguments()
    outdir = args.outdir
    prefix = args.prefix
    fasta_file = args.fasta
    bam_file = args.bam
    gap_len = args.gap_len
    keep_temp = args.keep_temp

    # Detect number of threads
    if args.threads is None:
        try:
            nproc = subprocess.run(
                "nproc", shell=True, check=True, text=True, capture_output=True
            )
            n_threads = int(nproc.stdout.strip())
        except subprocess.CalledProcessError:
            n_threads = os.cpu_count() or 1
    else:
        n_threads = args.threads
    print(f"{timenow()} Using {n_threads} CPUs for processing.")

    # If filter is requested, run BAM filter and exit
    if args.filter:
        chrA, chrB = args.filter
        os.makedirs(outdir, exist_ok=True)
        filter_bam(bam_file, chrA, chrB, outdir=outdir)
        return

    # Prepare output directories
    os.makedirs(outdir, exist_ok=True)
    combined_fasta_file = os.path.join(outdir, f"{prefix}_combined.fasta")
    bed_file = os.path.join(outdir, f"{prefix}.bed")
    genome_dir = os.path.join(outdir, "mini_genome")
    fastq_out = os.path.join(outdir, f"{prefix}.fastq")
    pri_bam = os.path.join(outdir, f"{prefix}_pri.bam")

    # Workflow
    read_fasta(fasta_file, combined_fasta_file, bed_file)
    check_star()
    mini_genome(combined_fasta_file, genome_dir, n_threads)
    bam_to_gapped_fastq(bam_file, fastq_out, gap_len=gap_len)
    map_to_minigenome(
        fastq_out,
        genome_dir,
        n_threads=n_threads,
        out_prefix=os.path.join(outdir, prefix),
    )

    if not keep_temp:
        cleanup(outdir, prefix)

    seglen_median = get_segment_lengths(pri_bam)
    crssant_birch(pri_bam, bed_file, outdir, n_threads, seglen_median)


if __name__ == "__main__":
    main()
