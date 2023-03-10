"""
seglendist_mod.py, Kongpan Li, 2022-04-15, modified from seglendist.py. Zhipeng Lu, 2020, zhipengluchina@gmail.com
This script calculates the size distribution of each segment in STAR alignments. 

Requires python3.

Example command for creating the list:
python seglendist.py align.sam align_seglen.pdf
"""

import matplotlib.pyplot as plt
import sys, re, matplotlib, matplotlib.pyplot as plt
from statistics import median


if len(sys.argv)< 3:
    print("Usage: python seglendist.py inputfile outputfile")
    sys.exit()

inputfile = open(sys.argv[1], 'r')
outputfile = sys.argv[2]


##part1: save the size list as space delimited numbers in one line in a file. 
sizelist = []
for line in inputfile:
    if line[0] == "@": continue
    CIGAR = line.split()[5]
    segs=[i.rstrip('0123456789') for i in CIGAR.split('N')]
    Mlens=[sum([int(i[:-1])for i in re.findall('\d+[M=X]',s)])for s in segs]
    sizelist += Mlens
inputfile.close()

print("Total number of segments", len(sizelist))
print("Segment length median", median(sizelist))

##part2: plot the arm size distribution
fig, ax = plt.subplots(figsize=(2.5,1.5))
plt.subplots_adjust(bottom=0.32)
plt.subplots_adjust(top=0.95)
plt.subplots_adjust(left=0.25)
plt.subplots_adjust(right=0.95)

n, bins, patches = plt.hist(sizelist, [x+0.5 for x in range(0,100)], \
                            histtype='step', cumulative=True, density=1) 
plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
plt.xlim(-10, 50)
plt.ylim(0, 1.1)
max_yticks = 5
yloc = plt.MaxNLocator(max_yticks)
ax.yaxis.set_major_locator(yloc)
ax.set_xlabel("segment length (nt)") #, fontsize=15
ax.set_ylabel("cumulative frequency") #, fontsize=15
plt.savefig(outputfile)
#plt.show()
           