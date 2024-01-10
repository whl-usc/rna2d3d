"""
contact:    wlee9829@gmail.com
date:       2024_01_05
python:     python3.10
script:     seqstats.py

This script is used to extract information
about the read counts, mapping statistics, 
and details of the gaptypes and gapfilter 
analyses from the CRSSANT pipeline.
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
    return str(datetime.now())[:-7]

def check_index(output_csv):
    try:
        df_output = pd.read_csv(output_csv)

    except FileNotFoundError:
        df_output = pd.DataFrame()

    # Defining the index column
    required_rows = ['Raw reads', 'Round 1', 'Input reads', 'Average input read length', 'Uniquely mapped reads', 'Multi-mapped reads', 'Too many mismatches', 'Too short reads', 'Other unmapped', 'Round 2', 'Input reads', 'Average input read length', 'Uniquely mapped reads', 'Multi-mapped reads', 'Too many mismatches', 'Too short reads', 'Other unmapped', 'Combined', 'Input reads', 'Average input read length', 'Uniquely mapped reads', 'Multi-mapped reads', 'Too many mismatches', 'Too short reads', 'Other unmapped', 'Gap type analysis', 'Total input alignment number', 'Contiuous alignments (no gaps)', 'Two-segment gapped alignments', 'Multi-segment gapped alignments', 'Other chimeric (different str/chr)', 'Overlapping chimeric (homotypic)', 'Bad alignments', 'Filtering for gap_1','Total single gapped alignments', 'Alignments with at least 1 good gap', 'Alignment with at least 2 good gaps', 'Filtering for gap_m', 'Total multiple gapped alignments', 'Alignments with at least 1 good gap', 'Alignment with at least 2 good gaps', 'Total number of gaps', 'Median gap length (all)', 'Median gap length (dist. plot)', 'Total number of segments', 'Median segment length']
    
    # Check for and make the index column if it does not already exist
    if 'Metrics' not in df_output.columns:
        df_output['Metrics'] = ''
    missing_rows = [row for row in required_rows if row not in df_output['Metrics'].values]
    df_output = pd.concat([df_output, pd.DataFrame({'Metrics': missing_rows})], ignore_index=True)

    # Save the DataFrame to output_csv
    return df_output.to_csv(output_csv, index=False)

def fastq_count(file_path):
    if os.path.exists(file_path):
        try:
            if file_path.endswith('.gz'):
                with gzip.open(file_path, "rt", encoding=None) as gz_file:
                    line_count_output = subprocess.check_output(['wc', '-l'], stdin=gz_file, universal_newlines=True)
            else:
                    line_count_output = subprocess.check_output(['wc', '-l'], stdin=file_path, universal_newlines=True)
            line_count = int(line_count_output.split()[0])
            read_count = line_count // 4
        except subprocess.CalledProcessError as e:
            read_count = 0
            print(timenow(), f" Error counting reads for {file_path}. Omiting count.")

    else:
        read_count = 0
        print(timenow(), f" Fastq file was not provided. Omitting count.")

    return read_count

def mapping_info(directory_name):
    mapping_data = pd.DataFrame()
    if os.path.exists(directory_name):
        for i in range(1, 3):
            info_lines = []
            log_file_path = os.path.join(directory_name, f'map_{i}', f'{directory_name}_{i}_Log.final.out')
            with open(log_file_path, 'r') as file:
                for line in file:
                    # Add the lines of data to include
                    if "Number of input reads" in line or \
                        "Average input read length" in line or \
                        "Uniquely mapped reads number" in line or \
                        "Number of reads mapped to multiple loci" in line or \
                        "Number of reads unmapped: too many mismatches" in line or \
                        "Number of reads unmapped: too short" in line or \
                        "Number of reads unmapped: other" in line:
                        info_lines.append(line.split()[-1])

            # Concatenate vertically
            df = pd.DataFrame(info_lines, columns=[f'Map_{i}'])
            mapping_data = pd.concat([mapping_data, df], axis=1)
        
        mapping_data[f'Map_1'] = pd.to_numeric(mapping_data[f'Map_1']).astype(int)
        mapping_data[f'Map_2'] = pd.to_numeric(mapping_data[f'Map_2']).astype(int)
        mapping_data['Combined'] = (mapping_data[f'Map_1'] + mapping_data[f'Map_2']).astype(int)
        mapping_data = pd.concat([pd.DataFrame(index=range(2)), mapping_data[f'Map_1'], pd.DataFrame(index=range(1)), mapping_data[f'Map_2'], pd.DataFrame(index=range(1)), mapping_data[f'Combined']]).reset_index()
        mapping_data = mapping_data.drop(columns=['index'])
        mapping_data.rename(columns={0: directory_name}, inplace=True)

    else:
        print(timenow(),f" Error in parsing mapping data. Check to see if there is a mistake in the directory_name and if Log.final.out files exist.")

    return mapping_data

def split_slurm(file_path, directory_name):
    if os.path.exists(file_path):
        # Break the slurm file into the appropriate sections
        with open(file_path, 'r') as file:
            section = []; current_section = []; inside_section = False
            start_pattern = "..... started STAR run" 
            end_pattern = "SA tag removal completed successfully."
            for line in file:
                if start_pattern in line:
                    current_section = [line]
                    inside_section = True
                elif inside_section:
                    current_section.append(line)
                    if end_pattern in line:
                        section.append(current_section)
                        inside_section = False

        directory_names = [re.search(r"name='([^']+)'", lists[74].strip()).group(1).replace("_pri_crssant.sam", "") for lists in section if re.search(r"name='([^']+)'", lists[74].strip())]
        slurm = pd.DataFrame()
        for group, section in zip(directory_names, section):
            if group == directory_name:         
                select_rows = list(range(29, 36)) + list(range(47, 50)) + list(range(60, 63)) + list(range(67, 72))
                slurm_select = pd.DataFrame(section).iloc[select_rows].applymap(lambda x: x.strip().split()[-1] if isinstance(x, str) else x)
                slurm_select.fillna('', inplace=True)

                ins_newline = [28, 36, 50]
                for row in ins_newline:
                    slurm_select.loc[row, :] = ''
                slurm = slurm_select.sort_index().reset_index(drop=True)
                slurm.rename(columns={0: directory_name}, inplace=True)
    else:
        slurm = pd.DataFrame(0, index=range(23), columns=[directory_name])
        print(timenow(),f" slurm.out file was not provided. Omitting count.")
        
    return slurm

def combine_outputs(directory_name, fastq_path, slurm_path, output_csv):
    check_index(output_csv)
    existing_data = pd.read_csv(output_csv)
    read_count = fastq_count(fastq_path)
    mapping_data = mapping_info(directory_name)
    slurm_data = split_slurm(slurm_path, directory_name)
    combined_df = pd.concat([mapping_data, slurm_data], axis=0, ignore_index=True)
    combined_df.loc[0, directory_name] = read_count

    if directory_name in existing_data.columns:
        print(timenow(),f" Data already exists for {directory_name}, updating values instead.")
        existing_data.loc[:, directory_name] = combined_df.loc[:, directory_name]
        existing_data.to_csv(output_csv, index=False)
    else:
        df_output = pd.concat([existing_data, combined_df], axis=1, ignore_index=False)
        df_output.to_csv(output_csv, index=False)
    print(timenow(),f" Job completed.")

def main():
    parser = argparse.ArgumentParser(
        prog="seqstats.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\

        This script is used to extract information from various files of the CRSSANT pipeline: 

        1. Number of reads from the pre-mapping fastq file.
        2. Output from two rounds of mapping (Log.final.out).
        3. Output from the gap_types and gap_filter analyses.

        NOTE: Arguments should be provided in the following order:

        1. Directory name
        2. Path to fastq file
        3. Path to slurm-*.out
        4. Output name for CSV
        """),
        usage="\npython3 %(prog)s [-h] [-m] directory_name fastq_path slurm_path output_csv")

    parser.add_argument('-m', '--metrics', action='store_true', help="""displays the definition of the metrics used in the script.""")
    
    args, remaining_args = parser.parse_known_args()
    
    if args.metrics:
        print("""
These output files have been previously described:

Integrated analysis of crosslink-ligation data for the detection of RNA 2D/3D conformations in vivo.
Lee WH, Zhang M and Lu Z. Methods in Molecular Biology (2023). 

Classification and clustering of RNA crosslink-ligation data reveal complex structures and homodimers.
Zhang M, Hwang IT, Li K, Bai J, Chen JF, Weissman T, Zou JY and Lu Z. Genome Research, 32: 1-18 (2022).

Other publications from the Zhipeng Lu Lab at USC have implemented such metrics in the corresponding analyses.

###########################################################################
###########################################################################

A brief description of the crucial metrics is supplied here.

##############################

Raw reads: Number of sequencing reads from the fastq file.

##############################

"Round 1": Generated output files after the first round of STAR mapping of the pre-processed fastq file. 
"Round 2": Generated output files after performing softreverse.py on the cont.sam file from Round 1 of STAR mapping as the input for the second round of STAR mapping.  
"Combined": The summation of the information from the two rounds of STAR mapping.

Input reads: 
Average input read length:
Uniquely mapped reads:
Multi-mapped reads:
Too many mismatches:
Too short reads:
Other unmapped:

##############################

Gap type analysis:

Total input alignment number:
Contiuous alignments (no gaps):
Two-segment gapped alignments:
Multi-segment gapped alignments:
Other chimeric (different str/chr):
Overlapping chimeric (homotypic):
Bad alignments:

##############################

Filtering for gap_1:
Filtering for gap_m:

Total single gapped alignments:
Alignments with at least 1 good gap:
Alignment with at least 2 good gaps:

Total number of gaps:
Median gap length (all):
Median gap length (dist. plot):
Total number of segments:
Median segment length:

###########################################################################
###########################################################################
""")
    
    parser.add_argument('directory_name', help="""Directory NAME containing map_1 and map_2 directories with the Log.final.out files. Should be identical to 'Outprefix' from the mapping steps. Providing a directory PATH will result in the last component of the path being used as the directory NAME. 

    CAUTION: 'directory_name' should be unique for each set of data. Using identical names will overwrite any existing data.""")
    
    parser.add_argument('fastq_path', help="""Absolute PATH to the fastq (.fastq) file. Input can be compressed (.fastq.gz). Type 'none' if it does not exist.""")
    
    parser.add_argument('slurm_path', help="""slurm-*.out absolute PATH.The file that is generated after the gaptype and gapfilter analyses. Type 'none' if 
    it does not exist.""")
    
    parser.add_argument('output_csv', help="""NAME for the output file. If the CSV name is reused, results will append to the rightmost column or update any existing values with the same directory_name.""")

    args = parser.parse_args(remaining_args)
    args.directory_name = os.path.basename(os.path.normpath(args.directory_name))
    combine_outputs(args.directory_name, args.fastq_path, args.slurm_path, args.output_csv)

if __name__ == "__main__":
    main()