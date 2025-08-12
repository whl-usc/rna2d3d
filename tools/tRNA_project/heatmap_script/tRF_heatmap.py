# /usr/bin/env python3

"""
Contact:    wlee9829@gmail.com
Date:       2025_07_30
Python:     python3.10
Script:     tRF_heatmap.py

Construct heatmap of start/stop locations based on MINTbase information.
"""

################################################################################
# Define script version
__version__ = "2.1.1"

# Version notes
__update_notes__ = """
2.1.1
    -   Changed condition without -r and -s to start from -1 instead of min.

2.1.0
    -   Auto-filter to remove "MT" mitochondrial tRNAs.

2.0.0
    -   Added -a, --amino-acid: filter by (Amino acid and anticodon).
    -   Added -r, --range: dynamic or manual range setting.
    -   Added -s, --scale: optionally scale by number of datasets.
    -   Added sanity check for column spacing or filtering errors from Excel.

1.1.0
    -   Added color scaling relative to number of studies.
    -   Changed output file type to SVG.
    -   Added -t, --type to filter tRF by #Type column.
    
1.0.0
    -   Initial commit of script, setup logic and functions as needed.
"""
################################################################################
# Import packages

import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import seaborn as sns
import sys

################################################################################
# Define functions


def read_csv(input_file, tRF_type="all", filter_aa=None):
    """
    Opens the MINTbase CSV and generates filtered dataframe.

        Type
        License Plate (sequence derived)
        MINTbase Alternative IDs (GRCh37 assembly-derived)
        tRNA number
        Amino acid and anticodon
        Chromosome
        Chromosome strand
        Chromosome start position
        Chromosome end position
        Start position relative to start of mature tRNA
        End position relative to start of mature tRNA
        Fragment sequence
        # Distinct anticodons
        # Instances in true nuclear tRNAs
        # Instances in true MT tRNAs
        # Instances in tRNA lookalikes in nucleus
        Fragment length
        D-loop overlap?
        Anticodon-loop overlap?
        Anticodon-triplet overlap?
        T-loop overlap?
        Exclusively within tRNA genes?
        Expressed (# of datasets)?
        Maximum RPM

    Parameters:
        input_file (str):   Name of the input file (should be in CSV).
    """
    if input_file.endswith(".csv"):
        mintbase_df = pd.read_csv(input_file)
        # print(mintbase_df.head())

    filtered_df = mintbase_df[
        [
            "#Type",
            "Fragment sequence",
            "Amino acid and anticodon",
            "Chromosome",
            "Start position relative to start of mature tRNA",
            "End position relative to start of mature tRNA",
            "Expressed (# of datasets)?",
        ]
    ].copy()

    filtered_df = filtered_df[filtered_df["Chromosome"] != "MT"]

    filtered_df.columns = [
        "class",
        "sequence",
        "amino_acid",
        "chromosome",
        "start_pos",
        "end_pos",
        "scale",
    ]

    # Map CLI type to actual type labels
    type_map = {
        "5half": "5'-half",
        "3half": "3'-half",
        "5prime": "5'-tRF",
        "3prime": "3'-tRF",
        "iTRF": "i-tRF",
        "all": None,
    }

    if tRF_type not in type_map:
        raise ValueError(f"Invalid --type argument: {tRF_type}")

    target_class = type_map[tRF_type]
    if target_class:
        filtered_df = filtered_df[filtered_df["class"] == target_class]

    # Optional amino acid filtering
    if filter_aa:
        if isinstance(filter_aa, str):
            filter_aa = [filter_aa]

        # Flatten all comma-separated values into a list
        raw_filters = []
        for aa_entry in filter_aa:
            raw_filters.extend(aa_entry.split(","))

        processed_filters = []
        for aa in raw_filters:
            aa = aa.strip()
            if len(aa) >= 3:
                formatted = (
                    aa[0].upper()
                    + aa[1:3].lower()
                    + (aa[3:].upper() if len(aa) > 3 else "")
                )
                processed_filters.append(formatted)

        def aa_filter(val):
            return any(
                val.startswith(filt) if len(filt) == 3 else val == filt
                for filt in processed_filters
            )

        filtered_df = filtered_df[filtered_df["amino_acid"].apply(aa_filter)]

    # Check for scale values "(Expressed (# of datasets)"
    if pd.api.types.is_string_dtype(filtered_df["scale"]):
        filtered_df["scale"] = filtered_df["scale"].str.extract(r"(\d+)")
    filtered_df["scale"] = pd.to_numeric(filtered_df["scale"], errors="coerce")
    if filtered_df["scale"].isna().any():
        print(
            f"Error: Non-numeric values found in 'scale' column after"
            f"conversion. Check {input_file} column headers..."
        )
        sys.exit(1)

    return filtered_df


def plot_heatmap(
    input_df, output_file=None, tRF_type="all", scaled=None, plot_range=None
):
    """
    Plots a heatmap of start vs end positions based on tRNA fragments.

    Parameters:
        input_df (pd.DataFrame): Includes "start_pos" and "end_pos" columns.
        output_file (str): If provided, saves the heatmap to this file (PNG).
    """
    # Ensure positions are integers
    df = input_df.copy()
    df["start_pos"] = pd.to_numeric(df["start_pos"], errors="coerce")
    df["end_pos"] = pd.to_numeric(df["end_pos"], errors="coerce")
    df["scale"] = pd.to_numeric(df["scale"], errors="coerce")

    df = df.dropna(subset=["start_pos", "end_pos", "scale"])
    df["start_pos"] = df["start_pos"].astype(int)
    df["end_pos"] = df["end_pos"].astype(int)

    # Create a pivot table (2D matrix of counts)
    if scaled:
        # Treat each row contribution as 1x scale, sum across matching positions.
        heatmap_data = (
            df.groupby(["end_pos", "start_pos"])["scale"].sum().unstack(fill_value=0)
        )
    else:
        heatmap_data = df.groupby(["end_pos", "start_pos"]).size().unstack(fill_value=0)

    if plot_range is not None:
        try:
            start_min, start_max, end_min, end_max = map(int, plot_range.split(","))
        except ValueError:
            print(
                f"[ERROR] Invalid format for --range. Use: start_min,start_max,end_min,end_max"
            )
            sys.exit(1)

    else:
        if not plot_range and scaled:
            start_min, start_max = -1, df["start_pos"].max()
            end_min, end_max = -1, df["end_pos"].max()
        else:
            # Auto-detect the range
            start_min, start_max = df["start_pos"].min(), df["start_pos"].max()
            end_min, end_max = df["end_pos"].min(), df["end_pos"].max()

    # Expand range slightly
    start_range = range(start_min, start_max + 1)
    end_range = range(end_min, end_max + 1)
    heatmap_data = heatmap_data.reindex(
        index=end_range, columns=start_range, fill_value=0
    )

    # Plot the data
    plt.figure(figsize=(8, 6))
    sns.set_style("white", {"axes.grid": False})
    ax = sns.heatmap(
        heatmap_data,
        cmap="Blues",
        cbar_kws={"shrink": 0.4},
        linewidths=0.10,
        linecolor="grey",
    )
    ax.invert_yaxis()

    # Add tick markers for the x and y
    ax.tick_params(
        axis="both",
        which="major",
        length=3,
        width=0.5,
        direction="out",
        color="black",
        labelsize=8,
        bottom=True,
        top=False,
        left=True,
        right=False,
    )

    # Add spines to axes
    for edge in ["left", "right", "top", "bottom"]:
        ax.spines[edge].set_visible(True)
        ax.spines[edge].set_color("black")
        ax.spines[edge].set_linewidth(0.75)

    type_map = {
        "5half": "5'-half",
        "3half": "3'-half",
        "5prime": "5'-tRF",
        "3prime": "3'-tRF",
        "iTRF": "i-tRF",
        "all": "All tRF",
    }

    label = type_map.get(tRF_type, tRF_type)

    # Label the plot, change the font sizes
    plt.title(f"{label} Start vs End Position Heatmap", fontsize=10, fontweight="bold")
    plt.xlabel("Start Position", fontsize=8)
    plt.ylabel("End Position", fontsize=8)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=480)
    else:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot tRF heatmaps from MINTbase CSV data.\n\n"
            "Example:\n"
            "tRF_heatmap.py MINTbase.csv -a=Gly -t=5half -r=1,72,15,87 -s"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("input_file", help="Input CSV file from MINTbase.")
    parser.add_argument(
        "-o", "--output", action="store_true", help="Write the output to a filename."
    )
    parser.add_argument(
        "-a",
        "--amino-acid",
        nargs="+",
        help="Amino acid filter by tRNA or codon (e.g., Gly or GlyCCC). "
        "One or more, comma-separated (e.g., Gly,Glu,ThrTGT))."
        "Case-insensitive.",
    )
    parser.add_argument(
        "-r",
        "--range",
        type=str,
        help="Manual setting of the range as comma-separated values: "
        "-r=start_min,start_max,end_min,end_max (e.g., -r=0,100,0,100).",
    )
    parser.add_argument(
        "-s",
        "--scaled",
        action="store_true",
        help="Scaling for plot by number of datasets per fragment.",
    )
    parser.add_argument(
        "-t",
        "--type",
        default="all",
        choices=["5half", "3half", "itrf", "5prime", "3prime", "all"],
        help="Filter by tRF class (default: all).",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    MINTbase_csv = args.input_file
    amino_acid = args.amino_acid
    plot_range = args.range
    scaled = args.scaled
    type = args.type

    filtered_df = read_csv(MINTbase_csv, type, amino_acid)

    if args.output:
        tag = "scaled" if scaled != "none" else "raw"
        output_filename = (
            os.path.splitext(MINTbase_csv)[0] + f"_{type}_{tag}_heatmap.svg"
        )
        plot_heatmap(filtered_df, output_filename, type, scaled, plot_range)
    else:
        plot_heatmap(filtered_df, None, type, scaled, plot_range)


if __name__ == "__main__":
    main()
