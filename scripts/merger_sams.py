"""
minjiezhang123@gmail.com    Mar,18,2020
python3

use to merger sam files from first and second round mapping
"""


#1. input and output setup
################################################################################
from datetime import datetime
def timenow(): return str(datetime.now())[:-7]
import sys
if len(sys.argv) < 3:
    print("Usage: python merger.py  first_sam  second_sam  label  samtype  outputfile")
    print("first_sam:   sam file from first round mapping")
    print("second_sam:  sam file from second round mapping")
    print("label:       label that will be added after readID, or 'none'")
    print("samtype:     gap1, gapm, homo, trans")
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
