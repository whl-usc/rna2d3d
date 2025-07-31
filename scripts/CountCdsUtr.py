"""
minjiezhang123@gmail.com    2020-01-16
Used to count the reads number of CDS, 5UTR and 3UTR region.
eg: python CountCdsUtr.py test.bam hg38mask14add_count.bed test

Only the longest transcript was used for annotation.
"""


#1. input and output setup
################################################################################
################################################################################
#this section sets up the input and output
import sys, numpy, os, re, itertools, random
from datetime import datetime
from multiprocessing import Process, Lock, Manager
from collections import Counter

if len(sys.argv) < 3:
    print("Usage:       python  CountCdsUtr.py inputbam  CdsUtr.bed  exclude_type  outputprefix")
    print("inputbam:    Bam file")
    print("CdsUtr.bed   Bed file of CDS and 5UTR, 3UTR regions")
    print("             eg: hg38mask14add_count.bed")
    print("exclude_type:eg: snoRNA,rRNA")
    print("outputprefix Outputprefix")
    print("")
    print("eg:  python CountCdsUtr.py   test.bam   hg38mask14add_count.bed   test")
    sys.exit()

inputbam = sys.argv[1]
annobed = sys.argv[2]
exclude_types = sys.argv[3].split(',')
outputprefix = sys.argv[4]
genes = {}; genesdict = {}; cdsdict = {}; intronsdict = {}; fiveUtrdict = {}; threeUtrdict = {}; biotypedict = {};  inputcount = 0 
################################################################################


#2. Processing bam file
##############################################################################
########################## Process bam to bed file ###########################
print(str(datetime.now())[:-7], "Process bam to bed file ...")
outputbed = outputprefix+'.bed'
os.system("bedtools bamtobed -split -i %s > %s" % (inputbam, outputbed))
os.system("awk '{if($0!~/^chr/){print \"chr\"$0} else {print $0}}' %s > %s" % (outputbed, 'temp.bed'))
os.system("sort -k1,1 -k2,2n -o %s %s" % (outputbed, 'temp.bed'))
#totalreads = len(open(outputbed,'rU').readlines())


## get total reads number
f= open(outputbed,'r')
readsdict = {}
for line in f:
    readID = line.split('\t')[3]
    readsdict[readID] = 1
f.close()
totalreads = len(readsdict)


######################## get length of CDS and UTR ###########################
print(str(datetime.now())[:-7], "Extract CDS/UTR length ...")
annofile = open(annobed,'r')
genelendict = {};   cdslendict = {};    intronlendict = {};    fiveUtrlendict = {};    threeUtrlendict = {}
for line in annofile:
    align = line.split('\t')
    length = int(align[2]) - int(align[1]) + 1
    genes[align[4]] = align[6].rstrip('\n')
    if align[4] not in genelendict: genelendict[align[4]] = 0
    if align[4] not in intronlendict: intronlendict[align[4]] = 0
    if align[4] not in cdslendict: cdslendict[align[4]] = 0
    if align[4] not in fiveUtrlendict: fiveUtrlendict[align[4]] = 0
    if align[4] not in threeUtrlendict: threeUtrlendict[align[4]] = 0
    genelendict[align[4]] += int(length)
    if align[5] == "CDS":   cdslendict[align[4]] += int(length)
    if align[5] == "Intron":   intronlendict[align[4]] += int(length)
    if align[5] == "five_utr":   fiveUtrlendict[align[4]] += int(length)
    if align[5] == "three_utr":   threeUtrlendict[align[4]] += int(length)
annofile.close()



##########################   Get overlap  ###########################
print(str(datetime.now())[:-7], "Get overlap between bed and anno file ...")
outputoverlap = outputprefix+'_overlap.bed'
os.system("/project/zhipengl_72/minjiez/software/bedtools2/bin/bedtools intersect -a %s -b %s -F 0.5 -wa -wb > %s" % (annobed, outputbed, outputoverlap))



#chr1    3216025 3216968 -       Xkr4    CDS     protein_coding  \
#chr1    3216327 3216348 A00208:91:HHKVTDRXX:1:2237:27543:4131   50      -
########################## Counting CDS, 5UTR, 3UTR reads ###########################
#readsdict = {}
print(str(datetime.now())[:-7], "Counting CDS, 5UTR, 3UTR reads ...")
inputbed = open(outputoverlap,'r')
for line in inputbed:
    align = line.split()
    gene, region, readID, biotype = align[4], align[5], align[10], align[6]
    inputcount += 1
    if not inputcount%1000000:
        print(str(datetime.now())[:-7], "Processed", inputcount, "reads ...")
    
    if biotype not in exclude_types:
        readsdict[readID] = 1
        if biotype not in biotypedict:  biotypedict[biotype]={}
        if gene not in genesdict: genesdict[gene]=0
        if gene not in cdsdict: cdsdict[gene]=0
        if gene not in intronsdict: intronsdict[gene]=0
        if gene not in fiveUtrdict: fiveUtrdict[gene]=0
        if gene not in threeUtrdict: threeUtrdict[gene]=0
    
        biotypedict[biotype][readID] = 1
        genesdict[gene] += 1
        if region == "CDS":
            cdsdict[gene]+=1
        if region == "Intron":
            intronsdict[gene]+=1
        if region == "five_utr":
            fiveUtrdict[gene]+=1
        if region == "three_utr":
            threeUtrdict[gene]+=1
inputbed.close()

#totalreads = len(readsdict)

output = open(outputprefix+'_RPKM.txt','w')
outstring="Gene" + '\t' + "Biotype" + '\t' + "gene_RPKM" + '\t' + "CDS_RPKM" + '\t' + "Intron_RPKM" + '\t' + "5UTR_RPKM" + '\t' + "3UTR_RPKM" + '\n'
output.write(outstring)
for i in sorted(genes):
    if i not in genesdict or str(genesdict[i]) == "0": generpkm = "0.01";
    else:
        generpkm = genesdict[i]/(genelendict[i]/1000 * int(totalreads)/1000000)
    
    if i not in cdsdict or str(cdsdict[i]) == "0": cdsrpkm = "0.01";
    else:
        cdsrpkm = cdsdict[i]/(cdslendict[i]/1000 * int(totalreads)/1000000)
    
    if i not in intronsdict or str(intronsdict[i]) == "0": intronrpkm = "0.01";
    else:
        intronrpkm = intronsdict[i]/(intronlendict[i]/1000 * int(totalreads)/1000000)
        
    if i not in fiveUtrdict or str(fiveUtrdict[i]) == "0":  fiveUtrrpkm = "0.01";
    else:
        fiveUtrrpkm = fiveUtrdict[i]/(fiveUtrlendict[i]/1000 * int(totalreads)/1000000)
    
    if i not in threeUtrdict or str(threeUtrdict[i]) == "0":    threeUtrrpkm = "0.01";
    else:
        threeUtrrpkm = threeUtrdict[i]/(threeUtrlendict[i]/1000 * int(totalreads)/1000000)
    
    outstring = i + '\t' +  genes[i]+ '\t' + str(generpkm) + '\t' + str(cdsrpkm) + '\t' + str(intronrpkm) + '\t' + str(fiveUtrrpkm) + '\t' + str(threeUtrrpkm) + '\n'
    output.write(outstring)
output.close()



output = open(outputprefix+'_RPM.txt','w')
outstring="Gene" + '\t' + "Biotype" + '\t' + "gene_RPM" + '\t' + "CDS_RPM" + '\t' + "Intron_RPM" + '\t' + "5UTR_RPM" + '\t' + "3UTR_RPM" + '\n'
output.write(outstring)

print('total reads:\t'+str(totalreads)+'\n')

for i in sorted(genes):
    if i not in genesdict or str(genesdict[i]) == "0": generpkm = "0.01";
    else:
        generpkm = genesdict[i]/int(totalreads)*1000000
    
    if i not in cdsdict or str(cdsdict[i]) == "0": cdsrpkm = "0.01";
    else:
        cdsrpkm = cdsdict[i]/int(totalreads)*1000000
    
    if i not in intronsdict or str(intronsdict[i]) == "0": intronrpkm = "0.01";
    else:
        intronrpkm = intronsdict[i]/int(totalreads)*1000000
        
    if i not in fiveUtrdict or str(fiveUtrdict[i]) == "0":  fiveUtrrpkm = "0.01";
    else:
        fiveUtrrpkm = fiveUtrdict[i]/int(totalreads)*1000000
    
    if i not in threeUtrdict or str(threeUtrdict[i]) == "0":    threeUtrrpkm = "0.01";
    else:
        threeUtrrpkm = threeUtrdict[i]/int(totalreads)*1000000
    
    outstring = i + '\t' +  genes[i]+ '\t' + str(generpkm) + '\t' + str(cdsrpkm) + '\t' + str(intronrpkm) + '\t' + str(fiveUtrrpkm) + '\t' + str(threeUtrrpkm) + '\n'
    output.write(outstring)
output.close()



output = open(outputprefix+'_ReadCount.txt','w')
outstring="Gene" + '\t' + "Biotype" + '\t' + "gene_Reads" + '\t' + "CDS_Reads" + '\t' + "Intron_Reads" + '\t' + "5UTR_Reads" + '\t' + "3UTR_Reads" + '\n'
output.write(outstring)
for i in sorted(genes):
    if i not in genesdict or str(genesdict[i]) == "0": generpkm = "0.01";
    else:
        generpkm = genesdict[i]
    
    if i not in cdsdict or str(cdsdict[i]) == "0": cdsrpkm = "0.01";
    else:
        cdsrpkm = cdsdict[i]
    
    if i not in intronsdict or str(intronsdict[i]) == "0": intronrpkm = "0.01";
    else:
        intronrpkm = intronsdict[i]
        
    if i not in fiveUtrdict or str(fiveUtrdict[i]) == "0":  fiveUtrrpkm = "0.01";
    else:
        fiveUtrrpkm = fiveUtrdict[i]
    
    if i not in threeUtrdict or str(threeUtrdict[i]) == "0":    threeUtrrpkm = "0.01";
    else:
        threeUtrrpkm = threeUtrdict[i]
    
    outstring = i + '\t' +  genes[i]+ '\t' + str(generpkm) + '\t' + str(cdsrpkm) + '\t' + str(intronrpkm) + '\t' + str(fiveUtrrpkm) + '\t' + str(threeUtrrpkm) + '\n'
    output.write(outstring)
output.close()

"""    
totalreads = 0    
for i in biotypedict:
    totalreads += len(biotypedict[i])
"""
    
output = open(outputprefix+'_biotype.txt','w')
for i in sorted (biotypedict):
    outstring = i + '\t' + str(len(biotypedict[i])) + '\t' + str(len(biotypedict[i])/totalreads) + '\n'
    output.write(outstring)
output.close()
 
 
os.system("rm %s %s %s" % (outputbed, 'temp.bed', outputoverlap))
