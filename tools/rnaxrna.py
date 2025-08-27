#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_08_21
Python:     python3.10
Script:     rnaxrna.py

Script re-processes pritrans gapped files to be mapped to a single mini-genome.
I) Check for STAR and samtools/pysam, input files as FASTA/BAM
Subfunction
1) Checks for overlapping regions "hotspots" for trans-chromosome reads.
2) Filters by user-defined chromosomes, based on hotspots.

II) Process input BAM and FASTA files for miniGenome generation
III) Re-map BAM to mini-genome
IV) Set up and run CRSSANT_birch to cluster newly mapped reads
"""

################################################################################
# Define version
__version__ = "2.0.0"

# Version notes
__update_notes__ = """
2.0.0
    -   Added filtering function for creating smaller BAM files.
    -   Added helper function to determine BAM overlaps.

1.0.0
    -   Set up logic and workflow.
    -   Set up a function to check for STAR aligner installation.
"""

################################################################################
# Import packages
from collections import defaultdict
from datetime import datetime
import argparse
import glob
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


def genome_overlaps(bam_path, window=15, min_reads=10):
    """
    Identify hotspots (clusters of reads) on each chromosome and 
    find cross-chromosomal hotspot pairs via read pairs.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    # 1. Collect read pairs
    read_positions = {}
    pairs = {}
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        qname = read.query_name
        chrom = bam.get_reference_name(read.reference_id)
        pos = read.reference_start
        if qname not in read_positions:
            read_positions[qname] = (chrom, pos)
        else:
            c1, p1 = read_positions[qname]
            c2, p2 = chrom, pos
            pairs[qname] = (c1, p1, c2, p2)

    bam.close()
    print(f"{timenow()} Collected {len(pairs)} read pairs...")

    # 2. Build hotspots per chromosome
    chrom_reads = defaultdict(list)
    for q, (c1, p1, c2, p2) in pairs.items():
        chrom_reads[c1].append((p1, q))
        chrom_reads[c2].append((p2, q))

    hotspots = {}
    print(f"{timenow()} Hotspots found per chromosome...\n")
    for chrom, positions in chrom_reads.items():
        positions.sort()
        clusters = []
        start, end, names = None, None, set()
        for pos, q in positions:
            if start is None:
                start, end, names = pos, pos, {q}
                continue
            if pos - end <= window:
                end = pos
                names.add(q)
            else:
                if len(names) >= min_reads:
                    clusters.append((start, end, names))
                start, end, names = pos, pos, {q}
        if len(names) >= min_reads:
            clusters.append((start, end, names))
        hotspots[chrom] = clusters
        if (len(clusters) > 0):
            print(f"\t{chrom}: {len(clusters)}")

    # 3. Link hotspots across chromosomes
    hotspot_pairs = defaultdict(set)
    for q, (c1, p1, c2, p2) in pairs.items():
        h1 = next(((c1, s, e) for s, e, ns in hotspots.get(c1, []) if q in ns), None)
        h2 = next(((c2, s, e) for s, e, ns in hotspots.get(c2, []) if q in ns), None)
        if h1 and h2 and h1 != h2:
            key = tuple(sorted([h1, h2]))
            hotspot_pairs[key].add(q)

    # 4. Print summary
    print(f"\nAlignment hotspot pairs (≥{min_reads} reads in each hotspot):\n")
    print(f"\tReads\tarm1_pos\tarm2_pos")

    if not hotspot_pairs:
        print("No cross-chromosomal hotspot pairs found")
    else:
        sorted_pairs = sorted(
            [(h1, h2, qnames) for (h1, h2), qnames in hotspot_pairs.items() if len(qnames) >= min_reads],
            key=lambda x: (x[0][0], x[0][1], x[1][0], x[1][1])
        )
        
        for h1, h2, qnames in sorted_pairs:
            c1, s1, e1 = h1
            c2, s2, e2 = h2
            print(f"\t{len(qnames)}\t{c1}:{s1}-{e1}\t{c2}:{s2}-{e2}")
    print("\n")


def filter_bam(input_bam, chrA, chrB, outdir="./"):
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

        print(f"{timenow()} Filtered BAM saved as: {os.path.basename(outname)}...\n")

    return outname


def read_fasta(file, out_fasta, out_bed, gap_len=10):
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
    combined_seq = ("N" * gap_len).join(genes[name] for name in gene_names)

    print(f"{timenow()} Preparing fasta genes: {combined_name}")
    with open(out_fasta, "w") as fasta:
        fasta.write(f">{combined_name}\n")
        for i in range(0, len(combined_seq), 60):
            fasta.write(combined_seq[i : i + 60] + "\n")
    try:
        subprocess.run(["samtools", "faidx", out_fasta], check=True)
        print(f"{timenow()} Indexed FASTA with samtools: {out_fasta}.fai")
    except subprocess.CalledProcessError:
        print(f"{timenow()} Failed to index FASTA with samtools: {out_fasta}")
        
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


def bam_to_gapped_fastq(bam_file, fastq_out):

    read_groups = defaultdict(list)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            read_groups[read.query_name].append(read)

    with open(fastq_out, "w") as out:
        for qname, reads in read_groups.items():
            # Identify left vs right arm by reference coordinate
            arms = {r.reference_start: r for r in reads if not r.is_unmapped}
            sorted_arms = sorted(arms.values(), key=lambda r: r.reference_start)

            if len(sorted_arms) == 1:
                # Only one arm mapped: pad with N equal to the other arm's
                r1 = sorted_arms[0]
                r2_seq = "N" * len(r1.query_sequence)
                r2_qual = "!" * len(r1.query_sequence)
            else:
                r1, r2 = sorted_arms[:2]
                r2_seq = r2.query_sequence
                r2_qual = "".join(chr(q + 33) for q in r2.query_qualities)

            seq = r1.query_sequence + r2.query_sequence
            qual = (
                "".join(chr(q + 33) for q in r1.query_qualities)
                + "".join(chr(q + 33) for q in r2.query_qualities)
            )

            out.write(f"@{qname}\n{seq}\n+\n{qual}\n")


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
        "--chimOutType", "WithinBAM", "HardClip"
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
    """
    Remove intermediate files generated by STAR and rnaxrna workflow.
    Handles FASTQ, log files, STAR BAM, and STAR temporary directories.
    Ignores NFS transient files gracefully.
    """
    files_to_remove = [
        f"{outdir}/{out_prefix}.fastq",
        f"./Log.out",
        f"{outdir}/{out_prefix}_Log.out",
        f"{outdir}/{out_prefix}_Log.progress.out",
        f"{outdir}/{out_prefix}_SJ.out.tab",
        f"{outdir}/{out_prefix}_Unmapped.out.mate1",
        f"{outdir}/{out_prefix}_Aligned.sortedByCoord.out.bam"
    ]

    for f in files_to_remove:
        try:
            if os.path.exists(f):
                os.remove(f)
        except FileNotFoundError:
            pass  # NFS or transient deletion issues

    # Remove STAR temporary directory
    star_tmp_dirs = [
        os.path.join(outdir, f"{out_prefix}__STARtmp")  # double underscore
    ]

    for tmp_dir in star_tmp_dirs:
        if os.path.exists(tmp_dir):
            try:
                shutil.rmtree(tmp_dir, ignore_errors=True)
            except Exception as e:
                print(f"Warning: failed to remove {tmp_dir}: {e}")


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
    print(f"{timenow()} Median seglen for CRSSANT_birch: {seglen_median}\n")


def parse_arguments():
    """
    Set up command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Process RNA-RNA interaction data for CRSSANT clustering.",
        usage=(" rnaxrna.py --bam BAM" 
            "\n\t[--overlap] MIN_READS"
            "\n\t[--filter CHR_A CHR_B] "
            "\n\t[--fasta FASTA [other flags]]"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Always required, positional argument
    parser.add_argument(
        "bam_file",
        help="Input BAM file containing trans-paired reads to convert into gapped FASTQ.",
    )

    # Optional function flags
    parser.add_argument(
        "--overlap",
        nargs="?",
        type=int,
        metavar="NUM_READS",
        const=10,
        default=None,
        help="SUB-FUNCTION 1: Checks BAM coverage per genomic region.",
    )
    parser.add_argument(
        "--filter",
        nargs=2,
        metavar=("CHR_A", "CHR_B"),
        help="SUB-FUNCTION 2: Filter BAM for reads connecting CHR_A and CHR_B.",
    )

    # Supplementary flags
    parser.add_argument(
        "--fasta",
        "-fa",
        help="Input FASTA file containing sequences to merge into a mini-genome.",
    )
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=None,
        help="Number of threads for STAR (auto-detected if not provided).",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        default="./",
        help="Output directory for mini-genome and mapping results.",
    )
    parser.add_argument(
        "--prefix", "-p", help="Prefix for STAR output."
    )
    parser.add_argument(
        "--gap-len",
        "-g",
        type=int,
        default=10,
        help="Number of Ns to insert between paired reads in recombined FASTQ.",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep intermediate files (default: remove).",
    )
    parser.add_argument(
        "--version", "-v", action="version", version=f"%(prog)s {__version__}"
    )
    args = parser.parse_args()

    # Conditional validation after parsing
    if not args.filter and not args.overlap and not args.fasta:
        parser.error("fasta is required when not using either of the subfunctions --filter/--overlap")

    return args
    

################################################################################
# Execute main
def main():
    args = parse_arguments()
    outdir = args.outdir
    fasta_file = args.fasta
    bam_file = args.bam_file
    if args.prefix:
        prefix = args.prefix
    else:
        prefix = os.path.basename(bam_file).replace("_filtered", "").rsplit(".", 1)[0]
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

    if args.overlap:
        genome_overlaps(bam_file, 100, args.overlap)
        return

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
    read_fasta(fasta_file, combined_fasta_file, bed_file, gap_len)
    check_star()
    mini_genome(combined_fasta_file, genome_dir, n_threads)
    bam_to_gapped_fastq(bam_file, fastq_out)
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