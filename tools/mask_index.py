"""
contact:    wlee9829@gmail.com
date:       2024_04_01
python:     python3.10
script:     mask_index.py

This script assesses indices based on a given sequence and input FASTA file. 
"""
# Define version
__version__="1.1.0"

# Verison notes
__update_notes__="""
1.1.0
    -   Added argument to define sequence orientations "identical | reverse |
        reverse_complement".
1.0.0
    -   Initial commit.

To do:
    -   Write output list of index locations for quicker reading.
    -   Change genome FASTA to be index based
    -   Add support for multithread processing. 
"""

# Import Packages
from collections import defaultdict
from datetime import datetime
import argparse
import os
import subprocess
import sys
import textwrap
import time

###########################################################################
# 1. Define common functions.
def timenow():
    """
    Current timestamp as a string.

    Returns:
        str: Current timestamp in format 'YYY-MM-DD HH:MM:SS'.
    """
    time = str(datetime.now())[:-7]

    return time

def rna2dna(input_file, orientation):
    """
    Converts the sequences in an input FASTA file into DNA (U to T).

    Returns:
        file: FASTA with the sequences converted to DNA. 
    """
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")
    # Extract file name without extension for output file.
    filename = input_file.split(".")[0]
    output_file = (f"{filename}.fa")

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.startswith(">"):
                if orientation == "identical":
                    modified_line = line.replace("U", "T")
                    outfile.write(modified_line)
                elif orientation == "reverse":
                    modified_line = line.replace("U", "T")
                    modified_line = modified_line[::-1]
                    outfile.write(modified_line)
                elif orientation == "reverse_complement":
                    complement = {'A': 'T', 'U': 'A', 'C': 'G', 'G': 'C'}    
                    dna_sequence = (''.join(complement.get(base, base) 
                       for base in reversed(line.strip())))
                    outfile.write(dna_sequence + "\n")
            else:
                outfile.write(line)

    return output_file

def read_fasta(filename):
    """
    Reads the sequences from a FASTA file.

    Returns:
        dict: {seq_name: sequence}
    """
    sequences = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                sequence_name = line.strip()[1:]
                sequences[sequence_name] = ''
            else:
                sequences[sequence_name] += line.strip()

    return sequences

def split_fasta(fasta_file):
    """
    Splits a FASTA file into separate files based on header lines using awk.

    Returns:
        files: FASTA with header sequences as names.
    """
    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(f"FASTA file '{fasta_file}' not found.")

    # Define the awk command to split the fasta file
    awk_command = ("awk '/^>/{f=substr($0,2)\".fasta\"; print > f; next} "
        "{print >> f}'")
    # Run the awk command using subprocess
    subprocess.run(f"cat {fasta_file} | {awk_command}", 
        shell=True, check=True)

    subprocess.run(f"mkdir split_fasta && mv *.fasta split_fasta", 
        shell=True, check=True)

def find_sequence_in_fasta(sequence, fasta_dir):
    """
    Checks the sequence against the fasta files in fasta_directory.

    Returns:
        List: filename | start | end
    """
    locations = {}
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fasta'):
            with open(os.path.join(fasta_dir, filename), 'r') as file:
                full_sequence = (''.join(line.strip() 
                    for line in file.readlines()[1:]))
                start = full_sequence.find(sequence)
                if start != -1:
                    end = start + len(sequence) - 1
                    locations[filename] = (start, end)
    return locations

###########################################################################
# 2. Main, define accepted arguments. 

def parse_args():
    """
    Parse the command-line arguments.
    """
    parser = argparse.ArgumentParser(
        prog="mask_index.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent("""\
###########################################################################
-L, --list = Specify the input list file containing sequence names.
-F, --fasta = Specify the input FASTA file containing genome sequences.
-O, --orientation = Specify the orientation for conversion:
                     - 'identical': Keep the sequences identical.
                     - 'reverse': Reverse the sequences.
                     - 'reverse_complement': Reverse complement the sequences.
                     Note: All RNA "U" are converted to "T".
###########################################################################
            """),
        usage="""\
Read a list of sequences (if RNA, converts to DNA with given orientation) and 
compare them to a genome fasta, writing out indices. 

\npython3 %(prog)s 
    -F=".fasta" 
    -L=".txt" 
    -O="identical|reverse|reverse_complement"
""")

    parser.add_argument("-F", "--fasta", 
        type=str, required=True, help="(default: .fasta)")
    parser.add_argument("-L", "--list", type=str, required=True,
        help="(default: .txt)")
    parser.add_argument("-O", "--orientation", 
        choices=["identical", "reverse", "reverse_complement"],
        default="identical", type=str, required=True,
        help="(default: identical)")    
    parser.add_argument('-V', '--version', action='version', 
        version=f'%(prog)s {__version__}\n{__update_notes__}', 
        help='Print version, update notes, exit.')
    ##########################################################################

    return parser.parse_args()

def main():

    # Split the genome FASTA
    if os.path.exists("split_fasta"):
        print(timenow(), "Directory 'split_fasta' already exists. "
            "Skipping the genome FASTA splitting step.")    
    else:
        print(timenow(), "Splitting the genome FASTA into separate files.")
        fasta = args.fasta
        split_fasta(fasta)

    print(timenow(), "Reading the sequences from the provided list "
        "and finding their locations relative to the genome FASTA.")    

    seq_list = args.list
    orientation = args.orientation
    if orientation == "identical":
        print(timenow(), "Converting all U to T and using sequences as is.")
    input_file = rna2dna(seq_list, orientation)
    sequences = read_fasta(input_file)
    current_directory = os.getcwd() 
    fasta_dir = os.path.join(current_directory, "split_fasta")

    for name, sequence in sequences.items():
        locations = find_sequence_in_fasta(sequence, fasta_dir)
        for filename, (start, end) in locations.items():
            if filename is not None:
                print(f"\n\tName:\t\t{name}")
                print(f"\tSequence:\t{sequence}")   
                print(f"\t{filename}:\tStart:\t{start}\tEnd:\t{end}")

    # Remove temporary file if necessary.
    filename = seq_list.split(".")[0]
    seqname_file = (f"{filename}.fa")
    os.remove(seqname_file)

if __name__ == "__main__":
    args = parse_args()
    main()