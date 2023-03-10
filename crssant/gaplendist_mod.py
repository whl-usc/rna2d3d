"""
gaplendist_mod.py Kongpan Li, 2022-10-02, likongpan@gmail.com, modified from gaplendist.py by Zhipeng Lu, 2020, zhipengluchina@gmail.com

gaplendist.py, a simple script to plot the length of gaps (N). The input file
can be any type of SAM. For example, from the direct output of STAR
we only have the distribution from normally gapped reads. 

the modified version 
1. directly creates the final distribution pdf without generating the gap list. if a gap list file is needed as source data, 
we can turn to use the original version, gaplendist.py
2. adds one more option of the xlim, x-axis limit, which can be, e.g., 100, 1000, ... 
we don't recommend using maximum gap length as the xlim for long genomes, like human genome, because the maxlen can be super long but with few reads, which
renders the gap length distribution very low resolution and informationless.
3. calculates the median of gap lengths.

Command example:
python gaplendist_mod.py input.sam input_gaplen.pdf xlim
"""

import sys, re, matplotlib
import matplotlib.pyplot as plt
from statistics import median

if len(sys.argv)< 3:
    print("Usage: python gaplendist.py inputfile outfile xlim")
    print("xlim: x-axis limit. if not provided, the default is 100")
    sys.exit()

if len(sys.argv) == 3:
    inputfile = open(sys.argv[1], 'r')
    xlim = 100

if len(sys.argv) == 4:
    inputfile = open(sys.argv[1], 'r')
    xlim = int(sys.argv[3])

sizelist = []
for line in inputfile:
    if line[0] == "@": continue
    align = line.split()
    cigar = align[5]
    sizelist += [int(i[:-1]) for i in re.findall('\d+N',cigar)]
inputfile.close()

plot_sizelist = []
for i in sizelist:
    if i <= xlim: plot_sizelist.append(i)

print("Total number of gaps", len(sizelist))
print("gap length median for all", median(sizelist))
print("gap length median of selected ones for the distribution plot", median(plot_sizelist))

fig, ax = plt.subplots(figsize=(2.5,1.5))
plt.subplots_adjust(bottom=0.32)
plt.subplots_adjust(top=0.95)
plt.subplots_adjust(left=0.25)
plt.subplots_adjust(right=0.95)
n, bins, patches = plt.hist(sizelist, [x+0.5 for x in range(0,xlim)], \
                            histtype='step', cumulative=True, density=1)
plt.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
plt.xlim(-10, xlim)
plt.ylim(0, 1.1)

max_yticks = 5
yloc = plt.MaxNLocator(max_yticks)
ax.yaxis.set_major_locator(yloc)
ax.set_xlabel("gap length (nt)") #, fontsize=15
ax.set_ylabel("frequency") #, fontsize=15
plt.savefig(sys.argv[2])
#plt.show()



