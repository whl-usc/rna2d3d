#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_12_18
Python:     python3.10
Script:     CountCdsUtr.py
"""
################################################################################
# Define version
__version__ = "1.3.0"

# Version notes
__update_notes__ = """
1.3.0
    -   Added automatic handling for paired-end reads. 
1.2.0
    -   Adjust biotypes to be by unique read counts.
1.1.0
    -   Rework counting logic to be stringent per base.
    -   Fix RPKM/RPM calculations to use all read counts.
1.0.0
    -   Refactor code from the original CountCdsUtr.py script.
"""
################################################################################

# Import Packages
from collections import defaultdict
from datetime import datetime
import argparse
import itertools
import os
import subprocess

################################################################################


# Define all sub-functions for processing
def get_time():
    return str(datetime.now())[:-7]


# 1. Process the BAM file
def bam2bed(input_bam, output_bed):
    """
    Convert BAM to BED (single-end or paired-end automatically)
    Returns number of unique fragments/reqds.
    """
    print(f"{get_time()}\tChecking for read pairs...")
    try: 
        check = subprocess.run(
            ["samtools", "view", "-c", "-f", "1", input_bam],
            capture_output=True, text=True, check=True
        )
        paired = int(check.stdout.strip()) > 0

    except subprocess.CalledProcessError:
        print(f"{get_time()}\tError checking BAM for paired-end reads. Assuming single-end.")
        paired = False

    temp_bam = input_bam
    if paired:
        print(f"{get_time()}\tPaired-end BAM detected. Name-sorting BAM for BEDPE...")
        # Create a temporary name-sorted BAM
        temp_bam = output_bed + ".namesorted.bam"
        subprocess.run(["samtools", "sort", "-n", "-o", temp_bam, input_bam], check=True)
        subprocess.run(["samtools", "index", temp_bam], check=True)

    temp_unsorted = output_bed + ".unsorted"

    print(f"{get_time()}\tConverting BAM to BED file using bedtools...")
    # Bedtools command
    bedtools_cmd = ["bedtools", "bamtobed", "-i", input_bam]
    if paired:

        bedtools_cmd.append("-bedpe")   # fragment-level BED for paired-end
        print(f"{get_time()}\tDetected paired-end BAM; using BEDPE output.")
    else:
        bedtools_cmd.append("-split")    # single-end BED
        print(f"{get_time()}\tDetected single-end BAM; using normal BED output.")

    # Run bedtools
    with open(temp_unsorted, "w") as out:
        subprocess.run(bedtools_cmd, stdout=out, check=True)

    # Sort BED
    sort_cmd = f"sort -k1,1 -k2,2n {temp_unsorted}"
    with open(output_bed, "w") as sorted_out:
        subprocess.run(sort_cmd, shell=True, stdout=sorted_out, check=True)

    read_ids = set()
    with open(output_bed) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 4:
                read_ids.add(fields[3])
    os.remove(temp_unsorted)
    print(f"{get_time()}\tProcessed {len(read_ids)} unique reads/fragments.")
    
    return len(read_ids)


# 2. Process the annotation file
def read_anno(annotation_bed):
    genelendict = defaultdict(int)
    cdslendict = defaultdict(int)
    intronlendict = defaultdict(int)
    fiveUtrlendict = defaultdict(int)
    threeUtrlendict = defaultdict(int)
    genes = {}
    biotypes = {}

    with open(annotation_bed, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            strand = fields[3]
            gene_id = fields[4]
            region_type = fields[5].lower()
            biotype = fields[6]

            length = end - start + 1
            genes[gene_id] = strand
            biotypes[gene_id] = biotype
            genelendict[gene_id] += length

            if region_type == "cds":
                cdslendict[gene_id] += length
            elif region_type == "intron":
                intronlendict[gene_id] += length
            elif region_type == "five_utr":
                fiveUtrlendict[gene_id] += length
            elif region_type == "three_utr":
                threeUtrlendict[gene_id] += length

    return {
        "genes": genes,
        "biotypes": biotypes,
        "gene_lengths": genelendict,
        "cds_lengths": cdslendict,
        "intron_lengths": intronlendict,
        "five_utr_lengths": fiveUtrlendict,
        "three_utr_lengths": threeUtrlendict,
    }


def get_overlap(annotation_bed, bed_file, outputprefix):
    output_overlap = f"{outputprefix}_overlap.bed"
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        annotation_bed,
        "-b",
        bed_file,
        "-F",
        "0.5",
        "-wa",
        "-wb",
        "-sorted",
    ]
    with open(output_overlap, "w") as fout:
        subprocess.run(cmd, stdout=fout, check=True)
    return output_overlap


def count_reads(overlap_file, exclude_types, gene_biotypes):
    counts = {
        "cds": defaultdict(int),
        "introns": defaultdict(int),
        "five_utr": defaultdict(int),
        "three_utr": defaultdict(int),
        "biotypes": {},
        "total_reads": 0,
        "gene_reads": defaultdict(int),
        "gene_read_sets": defaultdict(set),
    }

    exclude_set = set(exclude_types)

    with open(overlap_file) as f:
        for line in f:
            try:
                fields = line.strip().split("\t", 11)
                gene_id = fields[4]
                region_type = fields[5].lower()
                bio_type = fields[6]
                read_id = fields[10]
            except IndexError:
                continue

            if bio_type in exclude_types:
                continue

            counts["gene_reads"][gene_id] += 1
            counts["gene_read_sets"][gene_id].add(read_id)
            counts["biotypes"][gene_id] = bio_type
            counts["total_reads"] += 1

            if region_type == "cds":
                counts["cds"][gene_id] += 1
            elif region_type == "intron":
                counts["introns"][gene_id] += 1
            elif region_type == "five_utr":
                counts["five_utr"][gene_id] += 1
            elif region_type == "three_utr":
                counts["three_utr"][gene_id] += 1

    counts["gene_counts"] = dict(counts["gene_reads"])
    counts["gene_read_sets"] = counts["gene_read_sets"]

    return counts


def write_output(prefix, genes, biotypes, counts, total_reads, *dicts):
    rpkm_file = f"{prefix}_RPKM.txt"
    rpm_file = f"{prefix}_RPM.txt"
    count_file = f"{prefix}_ReadCount.txt"
    biotype_file = f"{prefix}_biotype.txt"

    genelen, cdslen, intronlen, fiveutrlen, threeutrlen = dicts[:5]
    gene_counts, cdsdict, introndict, fiveutrdict, threeutrdict = dicts[5:10]

    def calc_rpkm(count, length):
        return (
            0.01
            if count == 0 or length == 0
            else count / (length / 1000 * total_reads / 1_000_000)
        )

    def calc_rpm(count):
        return 0.01 if count == 0 else count / total_reads * 1_000_000

    with open(rpkm_file, "w") as f:
        f.write(
            "Gene\tBiotype\tgene_RPKM\tCDS_RPKM\tIntron_RPKM\t5UTR_RPKM\t3UTR_RPKM\n"
        )
        for gene in sorted(genes):
            biotype = biotypes.get(gene, "NA")
            line = (
                f"{gene}\t{biotype}\t"
                + "\t".join(
                    [
                        f"{calc_rpkm(gene_counts.get(gene, 0), genelen.get(gene, 0)):.2f}",
                        f"{calc_rpkm(cdsdict.get(gene, 0), cdslen.get(gene, 0)):.2f}",
                        f"{calc_rpkm(introndict.get(gene, 0), intronlen.get(gene, 0)):.2f}",
                        f"{calc_rpkm(fiveutrdict.get(gene, 0), fiveutrlen.get(gene, 0)):.2f}",
                        f"{calc_rpkm(threeutrdict.get(gene, 0), threeutrlen.get(gene, 0)):.2f}",
                    ]
                )
                + "\n"
            )
            f.write(line)

    with open(rpm_file, "w") as f:
        f.write("Gene\tBiotype\tgene_RPM\tCDS_RPM\tIntron_RPM\t5UTR_RPM\t3UTR_RPM\n")
        for gene in sorted(genes):
            biotype = biotypes.get(gene, "NA")
            line = (
                f"{gene}\t{biotype}\t"
                + "\t".join(
                    [
                        f"{calc_rpm(gene_counts.get(gene, 0)):.2f}",
                        f"{calc_rpm(cdsdict.get(gene, 0)):.2f}",
                        f"{calc_rpm(introndict.get(gene, 0)):.2f}",
                        f"{calc_rpm(fiveutrdict.get(gene, 0)):.2f}",
                        f"{calc_rpm(threeutrdict.get(gene, 0)):.2f}",
                    ]
                )
                + "\n"
            )
            f.write(line)

    with open(count_file, "w") as f:
        f.write(
            "Gene\tBiotype\tgene_Reads\tCDS_overlaps\tIntron_overlaps\t5UTR_overlaps\t3UTR_overlaps\n"
        )
        for gene in sorted(genes):
            biotype = biotypes.get(gene, "NA")
            line = (
                f"{gene}\t{biotype}\t"
                + "\t".join(
                    [
                        str(gene_counts.get(gene, 0)),
                        str(cdsdict.get(gene, 0)),
                        str(introndict.get(gene, 0)),
                        str(fiveutrdict.get(gene, 0)),
                        str(threeutrdict.get(gene, 0)),
                    ]
                )
                + "\n"
            )
            f.write(line)

    with open(biotype_file, "w") as f:
        f.write("Biotype\tgene_Reads\tFraction\tPercentage\n")

        read_biotype_map = defaultdict(set)
        for gene_id, reads in counts["gene_read_sets"].items():
            biotype = biotypes.get(gene_id, "unknown")
            for read_id in reads:
                read_biotype_map[read_id].add(biotype)

        biotype_counts = defaultdict(int)
        for biotypes in read_biotype_map.values():
            for biotype in biotypes:
                biotype_counts[biotype] += 1

        for biotype in sorted(biotype_counts):
            count = biotype_counts[biotype]
            fraction = count / total_reads if total_reads > 0 else 0
            f.write(f"{biotype}\t{count}\t{fraction:.6f}\t{fraction*100:.2f}%\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Counts reads in CDS, 5'UTR, and 3'UTR from BAM using BED and outputs files for biotype, ReadCount, RPKM, and RPM."
    )
    parser.add_argument("bam", help="Input BAM file (sorted and indexed)")
    parser.add_argument(
        "bed",
        help=(
            "Annotation BED file with regions labeled"
            "in column 4 (e.g., CDS, 5UTR, 3UTR)"
        ),
    )
    VALID_BIOTYPES = {
        "none",
        "lncRNA",
        "miRNA",
        "misc_RNA",
        "mt-tRNA",
        "protein_coding",
        "RNR",
        "rRNA",
        "scaRNA",
        "scRNA",
        "SNORA",
        "SNORD",
        "snRNA",
        "tRNA",
    }
    parser.add_argument(
        "exclude",
        type=lambda s: [x.strip() for x in s.split(",")],
        default=[],
        help="Exclude comma-separated bioType(s) (e.g., --exclude=miRNA,tRNA)",
    )
    parser.add_argument("outprefix", help="Output file prefix")

    return parser.parse_args()


def main():
    args = parse_args()
    bam_file = args.bam
    annotation_bed = args.bed
    exclude = set(args.exclude)
    outprefix = args.outprefix

    bed_file = outprefix + ".bed"
    number_reads = bam2bed(bam_file, bed_file)
    print(f"\t\t\t{number_reads} alignments processed.")
    anno_data = read_anno(annotation_bed)

    print(f"{get_time()}\tChecking for overlaps...")
    overlap_file = get_overlap(annotation_bed, bed_file, outprefix)
    counts = count_reads(overlap_file, exclude, anno_data["biotypes"])
    gene_counts = counts["gene_counts"]

    write_output(
        outprefix,
        set(gene_counts.keys()),
        counts["biotypes"],
        counts,
        number_reads,
        anno_data["gene_lengths"],
        anno_data["cds_lengths"],
        anno_data["intron_lengths"],
        anno_data["five_utr_lengths"],
        anno_data["three_utr_lengths"],
        gene_counts,
        counts["cds"],
        counts["introns"],
        counts["five_utr"],
        counts["three_utr"],
    )

    for temp_file in [bed_file, overlap_file]:
        if os.path.exists(temp_file):
            os.remove(temp_file)

    print(f"{get_time()}\tComplete. See files with prefix '{outprefix}'...")


if __name__ == "__main__":
    main()
