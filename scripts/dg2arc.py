'''
contact:    wlee9829@gmail.com
date:       2022_10_04
python:     python3.8
script:     dg2arc.py
    
This script is used to count base pairs of
duplex group (DGs) formed by CRSSANT, 
writing into a steploop arc file for IGV viewing.
'''

# Import packages
import sys, argparse, os, re, time, RNA
import numpy as np
from datetime import datetime

# Usage instructions
if len(sys.argv) < 4:
    print("Usage:       python dg_bp.py ref bedpe outname y/n DG_number")
    print("ref:         fasta file containing the sequence")
    print("bedpe:       bedpe file created after CRSSANT")
    print("outname:		output file name")
    print("y/n:			look at 'all' DG stem loops?")
    print("DG_number:	if 'n', which DG number?")
    sys.exit()

with open(sys.argv[1], 'r') as ref_fasta: 
    ref_tmp = ref_fasta.read().splitlines()
    ref = ref_tmp[1]
    ref_fasta.close()
with open(sys.argv[2], 'r') as bedpe:
    bedpe = bedpe.readlines()

outname = str(sys.argv[3])
which = str(sys.argv[4])

if which == "no":
	try:
		print(str(datetime.now())[:-7]+" Forming stemloops for DG "+str(sys.argv[5])+"\n")
	except:
		print(str(datetime.now())[:-7]+" No DG number specified. Exiting.")
		sys.exit()

####################################################################################################

def DNA_complement1(sequence):
	comp_dict = {"A":"T","T":"A","G":"C","C":"G","a":"t","t":"a","g":"c","c":"g"}
	sequence_list = list(sequence)
	sequence_list = [comp_dict[base] for base in sequence_list]
	string = ''.join(sequence_list)
	return string

def getfa(bed,ref):
	#bed: [chr1,12,25,gene,coverage,strand]
	RNAME,fastart,faend,STRAND = bed[0], bed[1], bed[2], bed[3]
	faseq = ''
	if STRAND == '+': faseq = ref[int(fastart)-1:int(faend)]
	if STRAND == '-': faseq = (DNA_complement1(ref[int(fastart)-1:int(faend)]))[::-1]
	return faseq

def RNAduplex(bedpe,ref):
		bedl = [align[0],int(align[1]),int(align[2]),align[8]]
		bedr = [align[3],int(align[4]),int(align[5]),align[9]]
		l_seq = getfa(bedl, ref);  r_seq = getfa(bedr, ref)

		#RNAduplex usage:
		#d.structure:    the basepair information
		#d.energy:       the energy information
		#d.i:            the last nucleotide of the first strand.
		#d.j:            the first nucleotide of the second strand. 

		fc = RNA.duplexfold(l_seq,r_seq)
		l_pair, r_pair = fc.structure.split('&')[0], fc.structure.split('&')[1] #basepair inf
		l_pairlen, r_pairlen = len(l_pair), len(r_pair) #basepair length
		l_seqlen, r_seqlen = len(l_seq), len(r_seq) #sequence length
		#print(fc.structure)
		#print(fc.i, fc.j)
		#print(l_seq, r_seq)
		#print(l_seqlen, r_seqlen)

		bedl_BP = ['','','']; bedr_BP = ['','','']
		bedl_BP[0] = bedl[0]; bedr_BP[0] = bedr[0];
		if bedl[3] == "+":
			bedl_BP[1] = bedl[1] + fc.i - l_pairlen
			bedl_BP[2] = bedl[1] + fc.i - 1
		else:
			bedl_BP[1] = bedl[1] + (l_seqlen - fc.i)
			bedl_BP[2] = bedl[1] + l_seqlen - (fc.i - l_pairlen + 1)

		if bedr[3] == "+":
			bedr_BP[1] = bedr[1] + fc.j - 1
			bedr_BP[2] = bedr[1] + fc.j + r_pairlen - 2
		else:
			bedr_BP[1] = bedr[1] + (r_seqlen - fc.j - r_pairlen) + 1
			bedr_BP[2] = bedr[1] + (r_seqlen - fc.j)
		
		return bedl_BP, bedr_BP

####################################################################################################

dg_list = []
for i in range(0,len(bedpe)):
	align = bedpe[i].strip('\n').split('\t')
	## US47    1       42      US47    4493    4530    A00208:91:HHKVTDRXX:2:2230:4960:33004   1       +       -
	dgs = RNAduplex(bedpe,ref)
	dg_list.append(dgs)

if sys.argv[4] == "y":
	with open(outname+'.bed', 'w') as outfile:
		outfile.write("track graphType=arc"+"\n")
		for i in range(0,len(dg_list)):
			DGi,start,end = i,dg_list[i][0][1],dg_list[i][1][2]
			steps = int(dg_list[i][0][2] - dg_list[i][0][1])
			for j in range(0,steps+1):
				stemloop_end = int(dg_list[i][0][1]+j*1),int(dg_list[i][1][2]-j*1)	
				outfile.write(str(dg_list[i][0][0])+'\t'+str(stemloop_end[0])+'\t'+str(stemloop_end[1])+'\t'+"stemloop"+str(i+1)+"\n")
	outfile.close()
	print(str(datetime.now())[:-7]+" Completed forming arcs for all DGs. Check for "+outname+".bed")

elif sys.argv[4] == "n":
	try:
		with open(outname+'.bed', 'w') as outfile:
			outfile.write("track graphType=arc"+"\n")
			i = int(sys.argv[5])-1
			DGi,start,end = i,dg_list[i][0][1],dg_list[i][1][2]
			steps = int(dg_list[i][0][2] - dg_list[i][0][1])
			for j in range(0,steps+1):
				stemloop_end = int(dg_list[i][0][1]+j*1),int(dg_list[i][1][2]-j*1)	
				outfile.write(str(dg_list[i][0][0])+'\t'+str(stemloop_end[0])+'\t'+str(stemloop_end[1])+'\t'+"stemloop"+str(i+1)+"\n")
		outfile.close()
		print(str(datetime.now())[:-7]+" Completed forming arcs for DG: "+sys.argv[5]+" Check for "+outname+".bed")
	except:
		print(str(datetime.now())[:-7]+" DGs not present in *.bedpe. Try a different number.")