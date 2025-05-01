#!/usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_02_19
Python:     python3.10
Script:     bedpetoarc.py

Update to the original bedpetobed12.py script, fixes error in conversion.
Converts bedpe records from CRSSANT to bed12 format, where (chrom1 == chrom2) 
and (strand1 == strand2). Adds RGB coloring based on DG scores.

bed12 (12 fields): 
chrom chromStart chromEnd name score strand thickStart thickEnd
itemRgb blockCount blockSizes blockStarts
"""

################################################################################
# Define version
__version__ = "2.0.0"

# Version notes
__update_notes__ = """
2.0.0
    -   Switched to sys.argv method for input files.
    -   Implemented color_arcs() to adjust RGB based on DG scores.
    -   Added optional --no_color flag to disable RGB arc coloring.
"""

################################################################################
# Import packages

import argparse
import math
import os

################################################################################
# Define sub-functions for processing

def bedpe_to_bed12(bedpe_path, add_color=True):
    """
    Converts BEDPE file to BED12 format with optional coloring.

    Args:
        bedpe_path (str): Path to input BEDPE file.
        add_color (bool): Apply RGB coloring if True.

    Returns:
        list: Converted BED12 records.
    """
    bed12_records = []

    with open(bedpe_path, 'r') as bedpefile:
        for line in bedpefile:
            record = line.strip().split()
            chrom = record[0]
            chromStart1, chromEnd1 = int(record[1]), int(record[2])
            chromStart2, chromEnd2 = int(record[4]), int(record[5])
            name = f"{record[6]}{record[7]}"
            score = record[8]
            strand = record[9]
            itemRgb = "0,0,0"

            blockSizes = str(
                f"{chromEnd1 - chromStart1},{chromEnd2 - chromStart2}")
            blockStarts = str(f"0,{chromStart2 - chromStart1}")

            bed12_record = [
                chrom, str(chromStart1), str(chromEnd2), 
                name, score, strand,
                str(chromStart1), str(chromStart1), 
                itemRgb, "2", blockSizes, blockStarts
            ]

            bed12_records.append(bed12_record)

    # Apply coloring unless explicitly disabled
    if add_color:
        bed12_records = color_arcs(bed12_records)

    print(
        f'Processed {len(bed12_records)} lines from '
        f'file: {(bedpe_path).split(".bedpe")[0]}.\n')

    return bed12_records


def color_arcs(bed12_records, max_intensity=200):
    """
    Assigns RGB colors to BED12 records based on DG score, sorting only within
    records that share the same base name (split from column 4).

    Args:
        bed12_records (list): List of BED12 records.
        max_intensity (int): Max RGB intensity (0-255) for highest DG score.

    Returns:
        list: BED12 records with updated RGB values.
    """
    from collections import defaultdict

    grouped_records = defaultdict(list)

    # Extract score from column 4 and group by base_name
    for rec in bed12_records:
        name_parts = rec[3].split(',')
        if len(name_parts) < 4:
            continue  # skip malformed entries
        base_name = ','.join(name_parts[:2])
        dg_score = float(name_parts[-1])
        rec[3] = ','.join(name_parts[:3])  # keep name + ID if needed
        rec[4] = str(dg_score)             # overwrite BED12 score column
        grouped_records[base_name].append(rec)

    # Define base colors to cycle through
    color_list = [
        (0, 0, 255),     # Blue
        (255, 0, 0),     # Red
        (0, 128, 0),     # Green
        (255, 0, 255)    # Purple
    ]

    updated_records = []
    gamma = 1.5

    for idx, (base_name, group) in enumerate(grouped_records.items()):
        base_color = color_list[idx % len(color_list)]
        scores = [float(rec[4]) for rec in group]
        max_log_score = math.log1p(max(scores)) if scores else 1

        for rec in group:
            raw_score = float(rec[4])
            log_score = math.log1p(raw_score)
            intensity_scale = (1 - (log_score / max_log_score) ** gamma)
            adj_color = tuple(int(c * intensity_scale) for c in base_color)
            rec[8] = f"{adj_color[0]},{adj_color[1]},{adj_color[2]}"
            updated_records.append(rec)

    return updated_records

def parse_args():
    """
    Parses command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments with BEDPE and BED12 file paths.
    """
    parser = argparse.ArgumentParser(
        description="Converts CRSSANT BEDPE file to BED12 format \
            for arc visualization, coloring enabled by default."
    )
    parser.add_argument("bedpe", help="Input BEDPE file summarizing DGs.")
    parser.add_argument("bed12", help="Output BED12 file with arcs.")
    parser.add_argument(
        "--no_color", action="store_false", dest="add_color",
        help="Disable RGB coloring based on DG scores."
    )

    return parser.parse_args()

def main():
    """
    Main function to execute the conversion script.
    """
    args = parse_args()

    if not os.path.exists(args.bedpe):
        raise FileNotFoundError(f"Input BEDPE file not found: {args.bedpe}")

    bed12_records = bedpe_to_bed12(args.bedpe, add_color=args.add_color)

    # Write to BED12 file
    outfile = args.bed12 if args.bed12.endswith(".bed") else f"{args.bed12}.bed"

    with open(outfile, 'w') as bed12file:
        bed12file.write(
            f'track name="{args.bed12}" graphType=arc '
            f'itemRgb="On" height=100\n'
        )
        for rec in bed12_records:
            bed12file.write("\t".join(rec) + "\n")

################################################################################
# Execute main

if __name__ == "__main__":
    main()