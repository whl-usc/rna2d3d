"""
contact:    wlee9829@gmail.com
date:       2024_02_05
python:     python3.10
script:     sam2igv.py

This script is used to convert a SAM file into a sorted BAM with index
for viewing using IGV.
"""

# Import Packages
import os
import subprocess
import sys

###########################################################################

def sam2igv(input_sam):
    # Generate output file names based on input prefix
    input_prefix = os.path.splitext(input_sam)[0]
    unsorted_bam = f"{input_prefix}_unsorted.bam"
    sorted_bam = f"{input_prefix}_sorted.bam"

    subprocess.run(['samtools', 'view', '-bS', '-o', unsorted_bam, input_sam])
    subprocess.run(['samtools', 'sort', '-o', sorted_bam, unsorted_bam])
    subprocess.run(['samtools', 'index', sorted_bam])
    subprocess.run(['rm', unsorted_bam])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sam2igv.py input.sam")
        sys.exit(1)

    input_sam = sys.argv[1]
    sam2igv(input_sam)
