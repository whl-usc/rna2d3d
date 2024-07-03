"""
contact:    wlee9829@gmail.com
date:       2024_01_05
python:     python3.10
script:     seqstats.py

This script is used to extract information about the read counts, mapping 
statistics, and details of the gaptypes and gapfilter analyses from the 
CRSSANT pipeline. Files should be organized according to the mapping shell 
script from the rna2d3d repository.  
"""

# Import Packages
from datetime import datetime
import argparse
import gzip
import os
import pandas as pd
import re
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

def check_index(output_csv):
    """
    Checks the existence of the specified CSV file and its structure. If the 
    file does not exist or is missing certain rows, it adds necessary rows to 
    maintain a consistent structure.

    Args:
        output_csv (str): The PATH to the CSV file.

    Returns:
        None
    """
    try:
        df_output = pd.read_csv(output_csv)

    except FileNotFoundError:
        df_output = pd.DataFrame()

    # Defining the index column
    required_rows = [
        'Raw reads',
        'Round 1',
        'Input reads',
        'Average input read length',
        'Uniquely mapped reads',
        'Multi-mapped reads',
        'Too many mismatches',
        'Too short reads',
        'Other unmapped',
        'Chimeric reads',
        'Round 2',
        'Input reads',
        'Average input read length',
        'Uniquely mapped reads',
        'Multi-mapped reads',
        'Too many mismatches',
        'Too short reads',
        'Other unmapped',
        'Chimeric reads',
        'Combined',
        'Input reads',
        'Average input read length',
        'Uniquely mapped reads',
        'Multi-mapped reads',
        'Too many mismatches', 
        'Too short reads',
        'Other unmapped',
        'Chimeric reads',
        'Gap type analysis (Round 1)',
        'Total input alignment number',
        'Continuous alignments (no gaps)',
        'Two-segment gapped alignments',
        'Multi-segment gapped alignments',
        'Other chimeric (different str/chr)',
        'Overlapping chimeric (homotypic)',
        'Bad alignments',
        'Gap type analysis (Round 2)',
        'Total input alignment number',
        'Continuous alignments (no gaps)',
        'Two-segment gapped alignments',
        'Multi-segment gapped alignments',
        'Other chimeric (different str/chr)',
        'Overlapping chimeric (homotypic)',
        'Bad alignments',
        'Filtering for gap_1',
        'Total single gapped alignments',
        'Alignments with at least 1 good gap',
        'Alignment with at least 2 good gaps',
        'Filtering for gap_m',
        'Total multiple gapped alignments',
        'Alignments with at least 1 good gap',
        'Alignment with at least 2 good gaps',
        'Total number of gaps',
        'Median gap length (all)',
        'Median gap length (dist. plot)',
        'Total number of segments',
        'Median segment length'
    ]
    
    # Check for and make the index column if it does not already exist
    if 'Metrics' not in df_output.columns:
        df_output['Metrics'] = ''
    missing_rows = [
        row for row in required_rows if row not in df_output['Metrics'].values]
    df_output = pd.concat(
        [df_output, pd.DataFrame({'Metrics': missing_rows})], ignore_index=True)

    # Save the DataFrame to output_csv
    return df_output.to_csv(output_csv, index=False)

def fastq_count(file_path):
    """
    Counts the number of reads in a FASTQ file.

    Args:
        file_path (str): Path to the FASTQ file.

    Returns:
        int: Number of reads in the FASTQ file.
    """
    # Try extracting fastq count data from either .fastq or .fastq.gz file. 
    # Divide total line numbers by 4 to get a count.
    if os.path.exists(file_path):
        try:
            with (gzip.open(file_path, "rt", encoding=None) if file_path.
                endswith('.gz') else open(file_path, 'r')) as file:
                    read_count = sum(1 for line in file if line.startswith('@'))
        except subprocess.CalledProcessError as e:
            read_count = 0
            print(timenow(), f" Error counting reads for {file_path}." 
                " Omiting count.")
    else:
        read_count = 0
        print(timenow(), f" Fastq file was not provided. Omitting count.")

    return read_count

def mapping_info(directory_name):
    """
    Extracts mapping information from Log.final.out files in the specified 
    directories. File should be in the map_* subfolders if using the rna2d3d 
    pipeline to perform analysis.

    Args:
        directory_name (str): Name of the directory containing map_1 and map_2.

    Returns:
        pd.DataFrame: DataFrame containing mapping information.
    """
    mapping_data = pd.DataFrame()
    if os.path.exists(directory_name):
        for i in range(1, 3):
            info_lines = []
            log_file_path = os.path.join(directory_name, f'map_{i}', 
                f'{directory_name}_{i}_Log.final.out')
            with open(log_file_path, 'r') as file:
                for line in file:
                    # Add the lines of data to include
                    if ("Number of input reads" in line or 
                        "Average input read length" in line or 
                        "Uniquely mapped reads number" in line or 
                        "Number of reads mapped to multiple loci" in line or 
                        "Number of reads unmapped: too many mismatches" in line 
                        or "Number of reads unmapped: too short" in line 
                        or "Number of reads unmapped: other" in line 
                        or "Number of chimeric reads" in line):
                            info_lines.append(line.split()[-1])

            # Concatenate vertically
            df = pd.DataFrame(info_lines, columns=[f'Map_{i}'])
            mapping_data = pd.concat([mapping_data, df], axis=1)
        
        mapping_data[f'Map_1'] = pd.to_numeric(
            mapping_data[f'Map_1']).astype(int)
        mapping_data[f'Map_2'] = pd.to_numeric(
            mapping_data[f'Map_2']).astype(int)
        mapping_data['Combined'] = (mapping_data[f'Map_1'] +
            mapping_data[f'Map_2']).astype(int)
        mapping_data.iloc[1, 2] /= 2

        mapping_data = pd.concat([pd.DataFrame(index=range(2)),
            mapping_data[f'Map_1'], pd.DataFrame(index=range(1)),
            mapping_data[f'Map_2'], pd.DataFrame(index=range(1)),
            mapping_data[f'Combined'], pd.DataFrame(index=range(1))])
        mapping_data.rename(columns={0: directory_name}, inplace=True)
        mapping_data = mapping_data.reset_index(drop=True)

    else:
        print(timenow(),f" Error in parsing mapping data. Check to see if "
            " there is a mistake in the directory_name, if you are in the "
            " directory where directory_name exists, and if Log.final.out "
            " files exist.")

    return mapping_data

def split_slurm(file_path, directory_name):
    """
    Reads a SLURM file and extracts relevant information for directory.

    Args:
        file_path (str): Path to the SLURM file.
        directory_name (str): Name of the directory to extract information.

    Returns:
        pd.DataFrame: DataFrame containing SLURM information for directory.
    """
    if not os.path.exists(file_path):
        print(timenow(), "slurm.out file was not provided. Omitting count.")
        return pd.DataFrame(0, index=range(21), columns=[directory_name])

    try:
        start_pattern = "Started gapanalysis.py ..." 
        end_pattern = "SA tag removal completed successfully."

        sections = []; current_section = []

        with open(file_path, 'r') as slurm_text:
            for line in slurm_text:
                current_section.append(line)
                if end_pattern in line:
                    sections.append(current_section)
                    current_section = []
        
        selected_section = None
        for section in sections:
            for line in section:
                if directory_name in line:
                    selected_section = section
                    break
            if selected_section is not None:
                break

        slurm_data = []
        for line in selected_section:
            # Add the lines of data to include
            if("Total input alignment number" in line or
                "Continuous alignments (no gaps)" in line or 
                "Two-segment gapped alignments" in line or 
                "Multi-segment gapped alignments" in line or 
                "Other chimeric (different str/chr)" in line or 
                "Overlapping chimeric (homotypic)" in line or 
                "Bad alignments" in line or 
                "Total input alignment number" in line or
                "Continuous alignments (no gaps)" in line or 
                "Two-segment gapped alignments" in line or 
                "Multi-segment gapped alignments" in line or 
                "Other chimeric (different str/chr)" in line or 
                "Overlapping chimeric (homotypic)" in line or 
                "Bad alignments" in line or 
                "  Total alignments" in line or 
                "Alignments with at least 1 good gaps" in line or 
                "Alignments with at least 2 good gaps" in line or 
                "Total number of gaps" in line or 
                "gap length median for all" in line or 
                "gap length median of selected ones for the distribution plot"
                 in line or "Total number of segments" in line or 
                "Segment length median" in line):
                    slurm_data.append(line.split()[-1])
        
        slurm = pd.DataFrame(slurm_data); slurm.columns = [directory_name]
        insert_blank_lines = [7, 15, 19]

        # Add blank lines to the DataFrame at specified indices
        for index in insert_blank_lines:
            slurm = pd.concat([slurm.iloc[:index], 
                pd.DataFrame([[' '] * len(slurm.columns)], 
                    columns=slurm.columns), 
                    slurm.iloc[index:]], ignore_index=True)

    except Exception as e:
        print(timenow(), f" {directory_name} was not found in the slurm.out file. Omitting count.")
        return pd.DataFrame(0, index=range(21), columns=[directory_name])
    
    return slurm

def combine_outputs(directory_name, fastq_path, slurm_path, output_csv):
    """
    Combines mapping, SLURM, and FASTQ count information, appends/updates CSV.

    Args:
        directory_name (str): Name of the directory.
        fastq_path (str): Path to the FASTQ file.
        slurm_path (str): Path to the SLURM file.
        output_csv (str): Path to the output CSV file.

    Returns:
        None
    """
    check_index(output_csv)
    existing_data = pd.read_csv(output_csv)
    read_count = fastq_count(fastq_path)
    mapping_data = mapping_info(directory_name)
    slurm_data = split_slurm(slurm_path, directory_name)
    combined_df = pd.concat(
        [mapping_data, slurm_data], axis=0, ignore_index=True)
    combined_df.loc[0, directory_name] = read_count
    
    if directory_name in existing_data.columns:
        print(timenow(),f" Data already exists for {directory_name}, "
            " updating values instead.")
        existing_data.loc[:, directory_name] = (
            combined_df.loc[:, directory_name])
        existing_data.to_csv(output_csv, index=False)
    else:
        df_output = pd.concat(
            [existing_data, combined_df], axis=1)
        df_output.to_csv(output_csv, index=False)

    print(timenow(),f" Job completed.")

def main():

    parser = argparse.ArgumentParser(
        prog="seqstats.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
###########################################################################
        
This script should be run in the parent directory of "directory_name". It will 
extract information from various files of the CRSSANT pipeline: 

-Number of reads from the pre-mapping fastq file.
-Output from two rounds of mapping (Log.final.out).
-Output from the gap_types and gap_filter analyses.

NOTE: Arguments should only be provided in the following order:

1. directory_name:  NAME of directory containing map_1 and directories with the 
                    Log.final.out files. Should be identical to 'Outprefix' 
                    from the mapping steps. Providing a directory PATH will 
                    result in the last component of the path being used as the 
                    directory_NAME. 

                    CAUTION: 'directory_name' should be unique for each 
                    set of data. Using identical names will overwrite 
                    existing data.

2. fastq_path:      Absolute PATH to the fastq (.fastq) file. Input can be 
                    compressed (e.g., fastq.gz). Type 'none' if does not exist.

3. slurm_path:      PATH of the file that is generated after the gaptype and 
                    gapfilter analyses. Type 'none' if it does not exist.

4. output_csv:      NAME for the output file. If the CSV name is reused, 
                    results will append to the rightmost column or update any 
                    existing values with the same directory_name.      

###########################################################################
"""),
        usage="\npython3 %(prog)s [-h] [-m] directory_name fastq_path "
            " slurm_path output_csv")

    parser.add_argument('-m', '--metrics', action='store_true', 
        help="""displays the definition of the metrics used in the script.""")
    
    args, remaining_args = parser.parse_known_args()
    
    if args.metrics:
        print("""

###########################################################################

These output files have been previously described:

Integrated analysis of crosslink-ligation data for the detection of RNA 2D/3D 
conformations in vivo. Lee WH, Zhang M and Lu Z. Methods in Molecular Biology 
(2023). 

Classification and clustering of RNA crosslink-ligation data reveal complex 
structures and homodimers. Zhang M, Hwang IT, Li K, Bai J, Chen JF, Weissman T, 
Zou JY and Lu Z. Genome Research, 32: 1-18 (2022).

Other publications from the Zhipeng Lu Lab at USC have implemented such metrics 
in their corresponding analyses.

###########################################################################

A brief description of the crucial metrics is supplied here.

Raw reads:  Number of sequencing reads from the fastq file.

"Round 1":  Output files after the first round of STAR mapping of the 
            pre-processed fastq file. 

"Round 2":  Output files after performing softreverse.py on the cont.sam file  
            from Round 1 of STAR mapping and using it as input for the second 
            round of mapping.  

"Combined": Summation of information from two rounds of STAR mapping.

Input reads:    Number of reads passed in to STAR mapper. 

Average input read length:  The average length of reads passed to STAR.

Uniquely mapped reads:  Number of reads mapped to a unique location in the 
                        reference genome.

Multi-mapped reads:     Reads mapped to multiple locations of the genome.

Too many mismatches:    Reads that have a number of mismatches exceeding the
                        threshold values specified in the STAR parameters. 

Too short reads:    Reads shorter than length specified in STAR parameters.

Other unmapped:     Unmappable reads to the reference genome.

Chimeric reads:     Reads aligned to different chromosomes.

###########################################################################

Gap type analysis:  Script used to differentiate reads into various "types". 
                    See the CRSSANT paper for detailed information.

Total input alignment number:   Total number of alignments (mapped reads) after 
                                STAR mapping.

Contiuous alignments (no gaps): Entire read mapped without gaps.

Two-segment gapped alignments:  Single gap between two segments.

Multi-segment gapped alignments: Multiple gaps between segments.

Other chimeric (different str/chr): Chimeric alignments involving different 
                                    chromosomes or strands.

Overlapping chimeric (homotypic):   Chimeric alignments within the same 
                                    chromosome or genomic structure, but with 
                                    overlapping segments.

Bad alignments: Alignments that are considered poor quality for various reasons.

###########################################################################

Filtering for gap_1/gap_m:  Applying a filter to remove splice alignments and 
                            short deletions (e.g., 1-2 nt gaps) 

Total single gapped alignments: Total number of alignments with a single gap.

Alignments with at least 1 good gap: Alignments w/ >= 1 high-quality gap.

Alignment with at least 2 good gaps: Alignments w/ >= 2 high-quality gaps. 

Total number of gaps: The overall count of gaps in the alignments.

Median gap length (all): The median length of all gaps in the alignments.

Median gap length (dist. plot): The median length of gaps used for distribution 
                                plotting (i.e., in the PDF output from the
                                mapping step)

Total number of segments:   Total count of segments in the alignments.

Median segment length:  Median length of all segments in the alignments.

###########################################################################
""")
    
    parser.add_argument('directory_name')
    parser.add_argument('fastq_path')
    parser.add_argument('slurm_path')
    parser.add_argument('output_csv')

    args = parser.parse_args(remaining_args)
    args.directory_name = os.path.basename(
        os.path.normpath(args.directory_name))
    combine_outputs(args.directory_name, 
        args.fastq_path, args.slurm_path, args.output_csv)

if __name__ == "__main__":
    main()
sys.exit()