'''
contact:    wlee9829@gmail.com
date:       2022_03_08
python:     python3.10
script:     sam2bamsort.py

This script is used to count convert SAM files
into BAM format, sorting, then indexing
for visualization in the IGV program.
'''

# Import packages
import sys, argparse, os, pysam, subprocess, time
from datetime import datetime
def timenow(): return str(datetime.now())[:-7]

# Usage instructions
if len(sys.argv) < 2:
    print("Usage:	python sam2bamsort.py sam")
    print("outname:	sam file name")
    sys.exit()

# File input and output names
sam=sys.argv[1]
prefix = str(os.path.splitext(sys.argv[1])[0])
view = str(prefix+".bam")
sort = str(prefix+"_sorted.bam")

# Run pysam
print(timenow()+" Starting samtools to view..."); pysam.view(sam, "-bS", "-o", view, catch_stdout=False) 
print(timenow()+" Starting samtools to sort..."); pysam.sort("-o", sort, view, catch_stdout=False)
print(timenow()+" Starting samtools to index..."); pysam.index(sort, catch_stdout=False)
print(timenow()+" Completed...Check "+view+" and "+sort+".bai in IGV.")