'''
contact:    wlee9829@gmail.com
date:       2023_11_28
python:     python3.10
script:     merger_sams.py
    
This script is used to merge SAM files
after the mapping and gaptypes.py analysis.
User defined range and RNA region.
'''
# Import Packages
from datetime import datetime
import sys

####################################################################################################

def timenow():
    return str(datetime.now())[:-7]

def process_sam_file(sam_file, label, sam_type, out_sam):
    for line in sam_file:
        if line.startswith("@"):
            out_sam.write(line)
        else:
            align = line.rstrip('\n').split()
            if len(align) > 0:
                if sam_type == "trans":
                    outsam.write(align[0]+'-'+label+'\t'+'\t'.join(align[1:])+'\n')
                else:
                    outsam.write(align[0]+'-'+label+'\t'+'\t'.join(align[1:19])+'\n')
# Usage instructions:
def main():
    if len(sys.argv) < 6:
        print("Usage: python merger_sams.py first_sam second_sam label sam_type outputfile")
        print("first_sam:   sam file from the first round mapping")
        print("second_sam:  sam file from the second round mapping")
        print("label:       label that will be added after readID, or 'none'")
        print("sam_type:    cont, gap1, gapm, homo, trans")
        print("outputfile:  the name of the output file")
        sys.exit()

    sam1 = open(sys.argv[1], "r")
    sam2 = open(sys.argv[2], "r")
    label = sys.argv[3]
    sam_type = sys.argv[4]
    out_sam = open(sys.argv[5], "w")

    print(timenow() + " Started merger_sams.py ...")

    for line in sam1:
        if line.startswith("@"):
            out_sam.write(line)
    process_sam_file(sam1, label, sam_type, out_sam)

    for line in sam2:
        if line.startswith("@"):
            out_sam.write(line)
    process_sam_file(sam2, label, sam_type, out_sam)

    sam1.close()
    sam2.close()
    out_sam.close()

if __name__ == "__main__":
    main()
    