'''
contact:    wlee9829@gmail.com
date:       2022_06_24
python:     python3.8
script:     dna2rna.py
    
This script is used to either convert 
between DNA and RNA sequences.

'''

# Import packages
import sys, time
from datetime import datetime
from Bio.Seq import Seq

# Usage instructions
if len(sys.argv) < 2:
    print("Usage:           python dna2rna.py sequence")
    print("sequence:        DNA sequence")
    sys.exit()

seq = Seq(sys.argv[1])
revseq = seq[::-1]
if "U" in seq:
    seq_type = "RNA"
else:
    seq_type = "DNA"

comp = seq.complement()
rcomp = seq.reverse_complement()

print("\n" + str(datetime.now())[:-7] + "\n")

print("           DNA or RNA: " + seq_type)
print("             Original: " + seq)
print("              Reverse: " + revseq)
if seq_type == "DNA":
    tran = seq.transcribe()
    print("                  RNA: " + tran)
else:
    tran = seq.back_transcribe()
    print("                  DNA: " + tran)

print("        Complementary: " + comp)
print("Reverse complementary: " + rcomp + "\n")