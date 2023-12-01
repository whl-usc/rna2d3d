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

# Usage instructions.
if len(sys.argv) < 5:
    print("Usage: python merger.py  sam_1  sam_2  label  samtype  outsam")
    print("sam_1:   sam file from first round mapping")
    print("sam_2:  sam file from second round mapping")
    print("label:       label that will be added after readID, or 'none'")
    print("samtype:     cont, gap1, gapm, homo, trans")
    print("outsam:     output file name")
    sys.exit()

sam1 = open(sys.argv[1],"r")
sam2 = open(sys.argv[2],"r")
label = sys.argv[3]
samtype = sys.argv[4]
outsam = open(sys.argv[5],"w")

#2. Processing the sam file
################################################################################
print(timenow()+" Started merger.py ...")
header='';
if label != 'none':
    for line in sam1:
        if line[0]=="@": outsam.write(line)
        else:
            align = line.rstrip('\n').split()
            if len(align)>0:
                if samtype=="trans":    outsam.write(align[0]+'-'+label+'\t'+'\t'.join(align[1:])+'\n')
                else:                   outsam.write(align[0]+'-'+label+'\t'+'\t'.join(align[1:19])+'\n')
    for line in sam2:
        if line[0]!="@":
            align = line.rstrip('\n').split()
            if len(align)>0:
                if samtype=="trans":    outsam.write(align[0]+'-'+label+'\t'+'\t'.join(align[1:])+'\n')
                else:                   outsam.write(align[0]+'-'+label+'\t'+'\t'.join(align[1:19])+'\n')

else:
    for line in sam1:
        if line[0]=="@": outsam.write(line)
        else:
            align = line.rstrip('\n').split()
            if len(align)>0:
                if samtype=="trans":    outsam.write(line)
                else:                   outsam.write('\t'.join(align[0:19])+'\n')
    for line in sam2:
        if line[0]!="@":
            align = line.rstrip('\n').split()
            if len(align)>0:
                if samtype=="trans":    outsam.write(line)
                else:                   outsam.write('\t'.join(align[0:19])+'\n')

sam1.close()
sam2.close()
outsam.close()