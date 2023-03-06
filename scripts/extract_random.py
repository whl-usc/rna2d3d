"""
This script can be used to extract certain number of target input file.

Usage:
python extract_random.py  type  input.sam  ReadsNum  output.sam
"""

#1. input and output setup
################################################################################
################################################################################
#this section sets up the input and output
import sys, os, re, random
from datetime import datetime
from itertools import combinations

if len(sys.argv) < 3:
    print('\n'+"Usage:  ")
    print("python  extract_random.py  inputfile  header_start_symbol  Number  output")
    print("")
    print("inputfile:           input file, can be sam or other type of files")
    print("header_start_symbol: eg: @, lines starting with @ will be skipped")
    print("Number:              random number")
    print("output:              output file name")
    sys.exit()

inputsam = open(sys.argv[1],'r')
symbol = str(sys.argv[2])
ReadsNum = int(sys.argv[3])
output = open(sys.argv[4],'w')
inputcount = 0;
################################################################################


#2 random select certain lines of input file
################################################################################
print(str(datetime.now())[:-7], "processing input file")
header = ''; lines=[]; 
for line in inputsam:
    if line[0]==symbol:
        header += line; continue
    inputcount += 1
    if not inputcount%100000:
        print(str(datetime.now())[:-7], "Processed", inputcount, "lines ...")
    
    lines.append(line)
inputsam.close()


# output 
output.write(header)
outputlines = random.sample(lines,ReadsNum)
output.write(''.join(outputlines))
output.close()
