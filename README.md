<div id="top"></div>

<!-- PROJECT NAME -->
# rna2d3d

The rna2d3d repository contains scripts used for analyzing data generated from RNA crosslinking, proximity-ligation, high-throughput sequencing experiments.

<!-- ABOUT THE PROJECT -->
## About the Project

Advancements in the field of RNA structure and interaction studies have enabled high throughput analysis of RNAs, revealing base pairing and higher order interactions. Despite the progress on the experimental methods for studying RNA secondary and intermolecular complexes, computational analysis of such data generated has remained a challenge, typically limited to specific studies. One subset of such experimental methods, the crosslinking- and proximity-ligation based approaches generates non-continuous reads that indicate specific RNA structures and interactions. Here, we detail a computational pipeline for analyzing the non-continuous "gapped" reads from a variety of methods that employ the crosslinking-ligation principle (*see* [PARIS](https://pubmed.ncbi.nlm.nih.gov/27180905/), [PARIS2](https://www.nature.com/articles/s41467-021-22552-y), [LIGR](https://pubmed.ncbi.nlm.nih.gov/27184080/), [SPLASH](https://pubmed.ncbi.nlm.nih.gov/27184079/), [SHARC](https://www.nature.com/articles/s41467-022-28602-3)). Briefly, the pipeline functions as follows. First, raw data are processed to remove sequencing adapters and PCR duplicates (and if applicable. separated by barcode). Next, all reads are permissively mapped to a custom reference genome using STAR. Then, primary alignments are extracted, classified, and assembled to duplex groups (DGs), non-overlapping (NGs), or tri-segment groups (TGs). Finally, DGs can be visualized using IGV, converted to stemloop groups, and have endpoints (read arm stopping points) determined. Furthermore, the nucleotide frequncy or SHAPE reactivity of DG arms can be determined. 

<!-- GETTING STARTED -->
# Getting Started

## Prerequisites

### **Hardware**

A high-performance compute (HPC) cluster with 64-bit Unix-based operating system, x86-64 compatible processer, and at least 32 GB of RAM is recommended for completing this protocol in its entirety\*. Disk space requirements are dependent on the size of raw sequencing reads (fastq) and subsequent aligned and processed read (sam/bam) files. Mapping and CRSSANT steps require an HPC cluster with > 30 GB memory. Post-mapping analysis steps can be performed locally on a local system >= 4 compute threads and 8 GB of RAM\*\*. Note that large datasets require more memory and upwards of one week for mapping and CRSSANT duplex group (DG) and non-overlapping group (NG) assembly. 

\* Operations described for pre-processing and mapping steps were performed on an HPC cluster (2.60 GHz 8-core Intel Xeon 2640v3, 32 GB, CentOS)

\*\* DG analysis steps on a standard desktop (3.3 GHz 4-core Intel Core i5, 8 GB, Ubuntu 20.04) 

### **Software**
1. [bedtools (v2.22+)](https://github.com/arq5x/bedtools2)
2. [CRSSANT](https://github.com/zhipenglu/CRSSANT/archive/refs/heads/master.zip)
3. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
4. [icSHAPE](https://github.com/qczhang/icSHAPE)
5. [IGV](https://software.broadinstitute.org/software/igv/download)
6. [Python3.x+](https://www.python.org/downloads/)
    * biopython
    * maptlotlib
    * NetworkX (v2.1+)
    * numpy
    * pandas
    * pysam
    * SciPy
    * seaborn
7. [samtools (v1.1+)](https://www.htslib.org/download/)
8. [STAR-2.7.1a](https://github.com/alexdobin/STAR/archive/refs/tags/2.7.1a.tar.gz)
9. [Trimmomatic-0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- PIPELINE STEPS -->
# Data pre-processing

Post-sequencing data is pre-processed to remove sequencing adapters, short reads, and PCR duplicates. Reads are separate by barcode, if applicable. Upon successful completion of these steps, sequencing data should be processed and suitable for mapping to a reference genome.

## STEP 1 : 3'-end adapter trimming

```java -jar trimmomatic-0.36.jar SE -threads 16 -phred33 x.fastq.gz x.trim3.fastq ILLUMINACLIP:P6SolexaRC35.fa:3:20:10 SLIDINGWINDOW:4:20 MINLEN:18 ```

## STEP 2 : Remove PCR duplicates

```perl readCollapse x_trim3.fastq x_trim3_nodup.fastq ```

## STEP 3 : Split by barcode

```perl splitFastq.pl -U x_trim3_nodup.fastq -l CGTGAT:R701 -b 6:6 -s ```

## STEP 4 : 5'-end adapter trimming

```java -jar trimmomatic-0.36.jar SE -threads 16 -phred33 bcnn.fastq x_trim_nodup_bcnn.fastq HEADCROP:17 MINLEN:20 ```

## STEP 5 : FastQC 

Visuallze inspect the quality of the pre-processed data using FastQC. The "per base sequence quality plot" average quality score should be high across all read positions (>=30). Adapter content should be close to 0.

<p align="right">(<a href="#top">back to top</a>)</p>

# STAR mapping to reference genome

After data pre-processing, reads must be mapped to a reference genome. Primary alignments can then be classified into different types and filtered for splice junction and short gapped reads. Reads can also be plotted for gap and segment lengths and basic read statistics determined. 

## STEP 1 : STAR mapping (permissive)

```
/STAR-2.7.1a/bin/Linux_x86_64/STAR \
--runThreadN 8 \
--runMode alignReads \
--genomeDir StaridxPath \
--readFilesIn Fastq \
--outFileNamePrefix Outprefix_1_ \
--genomeLoad NoSharedMemory \
--outReadsUnmapped Fastx \
--outFilterMultimapNmax 10 \
--outFilterScoreMinOverLread 0 \
--outSAMattributes All \
--outSAMtype BAM Unsorted SortedByCoordinate \
--alignIntronMin 1 \
--scoreGap 0 \
--scoreGapNoncan 0 \
--scoreGapGCAG 0 \
--scoreGapATAC 0 \
--scoreGenomicLengthLog2scale -1 \
--chimOutType WithinBAM HardClip \
--chimSegmentMin 5 \
--chimJunctionOverhangMin 5 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreDropMax 80 \
--chimNonchimScoreDropMin 20 \
```

## STEP 2 : Extract primary alignments

```
samtools view -h x_Aligned.sortedByCoord.out.bam | awk '$1~/^@/ || NF<21' > x_nonchimeric_temp.sam
samtools view -h x_Aligned.sortedByCoord.out.bam | awk '$1!~/^@/ && NF==21'> x_chimeric_temp.sam'
samtools view -bS -F 0x900 -o x_nonchimeric_pri.bam x_nonchimeric_temp.sam
samtools view -h x_nonchimeric_pri.bam > x_nonchimeric_pri.sam
cat x_nonchimeric_pri.sam x_chimeric_temp.sam > x_1_pri.sam
rm -f x_nonchimeric_temp.sam x_chimeric_temp.sam x_nonchimeric_pri.bam x_nonchimeric_pri.sam
```

## STEP 3 : Classify primary reads

```python3 gaptypes.py x_1_pri.sam x_1_pri -1 15 1 ```

## STEP 4 : Rearrange softclipped continuous reads

```python softreverse.py x_1_pricont.sam softrev.fastq ```

## STEP 5 : STAR mapping (second round, softrev.fastq as input)

```
/STAR-2.7.1a/bin/Linux_x86_64/STAR \
--runThreadN 8 \
--runMode alignReads \
--genomeDir StaridxPath \
--readFilesIn softrev.fastq \
--outFileNamePrefix Outprefix_1_ \
--genomeLoad NoSharedMemory \
--outReadsUnmapped Fastx \
--outFilterMultimapNmax 10 \
--outFilterScoreMinOverLread 0 \
--outSAMattributes All \
--outSAMtype BAM Unsorted SortedByCoordinate \
--alignIntronMin 1 \
--scoreGap 0 \
--scoreGapNoncan 0 \
--scoreGapGCAG 0 \
--scoreGapATAC 0 \
--scoreGenomicLengthLog2scale -1 \
--chimOutType WithinBAM HardClip \
--chimSegmentMin 5 \
--chimJunctionOverhangMin 5 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreDropMax 80 \
--chimNonchimScoreDropMin 20 \
```

## STEP 6 : Extract primary alignments 

```
samtools view -h x_2_Aligned.sortedByCoord.out.bam | awk '$1~/^@/ || NF<21' > x_2_nonchimeric_temp.sam
samtools view -h x_2_Aligned.sortedByCoord.out.bam | awk '$1!~/^@/ && NF==21'> x_2_chimeric_temp.sam'
samtools view -bS -F 0x900 -o x_2_nonchimeric_pri.bam x_2_nonchimeric_temp.sam
samtools view -h x_2_nonchimeric_pri.bam > x_2_nonchimeric_pri.sam
cat x_2_nonchimeric_pri.sam x_2_chimeric_temp.sam > x_2_pri.sam
rm -f x_2_nonchimeric_temp.sam x_2_chimeric_temp.sam x_2_nonchimeric_pri.bam x_2_nonchimeric_pri.sam
```

## STEP 7 : Classify primary reads

```python3 gaptypes.py x_2_pri.sam x_2_pri -1 15 1 ```

## STEP 8 : Combine output from both rounds of STAR mapping

```
python3 merger_sams.py x_1_prigap1.sam x_2_prigap1.sam gap1 x_prigap1.tmp
python3 merger_sams.py x_1_prigapm.sam x_2_prigapm.sam gapm x_prigapm.tmp
python3 merger_sams.py x_1_prihomo.sam x_2_prihomo.sam homo x_prihomo.tmp
python3 merger_sams.py x_1_pritrans.sam x_2_pritrans.sam trans x_pritrans.tmp

samtools view -H x_1_Aligned.sortedByCoord.out.bam > header

cat header x_prigap1.tmp > x_prigap1.sam 
cat header x_prigapm.tmp > x_prigapm.sam 
cat header x_prihomo.tmp > x_prihomo.sam 
cat header x_pritrans.tmp > x_pritrans.sam
rm -f header *.tmp
```

## STEP 9 : Filter splices and short gaps

```
python gapfilter.py Gtf x_prigap1.sam x_prigap1_filtered.sam 13 yes
python gapfilter.py Gtf x_prigapm.sam x_prigapm_filtered.sam 13 yes
```

## STEP 10a : Primary read statistics (Gap length)

```
cat x_prigap1_filtered.sam x_prigapm_filtered.sam > x_prigaps_filtered.sam
python gaplendist.py x_prigaps_filtered.sam sam x.list all
python gaplendist.py x.list list x_gaplen.pdf all
rm -f x_prigaps_filtered.sam x.list
```

## STEP 10b : Primary read statistics (Arm length)

```
cat x_prigap1_filtered.sam x_prigapm_filtered.sam x_pritrans.sam > x_prifiltered.sam
python seglendist.py x_prifiltered.sam sam x_prifiltered.list
python seglendist.py x_prifiltered.list list x_seglen.pdf
rm -f x_prifiltered.sam x_prifiltered.list
```

## STEP 11 : Primary read statistics (Abundance of reads mapped to genes)

```
samtools view -bS -o x_prigap1.bam x_prigap1.sam
python CountCdsUtr.py x_prigap1.bam CdsUtr.bed none x_pri
```

<p align="right">(<a href="#top">back to top</a>)</p>

# CRSSANT assembly

CRSSANT optimizes short-read mapping and clusters gap1 and trans alignments into duplex groups (DGs) and non-overlapping groups (NGs). 

## STEP 1 : Combined filtered gap1 and trans SAM files

```python merger.py x_prigap1_filtered.sam x_pritrans.sam x_pri_crssant.sam ```

## OPTIONAL : Tag reads with unique identifier

```
awk '$0!~/^@/ {printf $1"-TAG1\t"; for(i=2;i<=NF;++i) printf $i "\t"; printf"\n"}' x_pri_crssant.sam > x_pri_crssant_TAG1_tmp.sam
awk '$0!~/^@/ {printf $1"-TAG2\t"; for(i=2;i<=NF;++i) printf $i "\t"; printf"\n"}' x_pri_crssant.sam > x_pri_crssant_TAG2_tmp.sam
samtools view -H x_pri_crssant.sam > header
cat header x_pri_crssant_TAG1_tmp.sam x_pri_crssant_TAG2_tmp.sam > x_pri_crssant.sam
rm -f header x_pri_crssant_TAG1_tmp.sam x_pri_crssant_TAG2_tmp.sam
```

## STEP 2 : Build "genesfile" in BED format

This step can be performed manually by following the tab delimited format as described in the CRSSANT manual or paper. Alternatively, if the CountCdsUtr.py script is used, a "ReadCount.txt" file will be generated and can be input into genes2bed.py.

```python genes2bed.py x_ReadCount.txt CdsUtr.bed min_reads outname```

## STEP 3 : Generate bedgraphs

```
bedtools genomecov -bg -split -strand + -ibam x_pri_crssant.bam -g staridxPath/chrNameLength.txt > x_plus.bedgraph
bedtools genomecov -bg -split -strand - -ibam x_pri_crssant.bam -g staridxPath/chrNameLength.txt > x_minus.bedgraph
```

## STEP 4 : Intialize CRSSANT for DG and NG assembly

```python crssant.py -cluster cliques -t_o 0.2 -out ./ x_pri_crssant.sam $Bed x_plus.bedgraph,x_minus.bedgraph```

## STEP 5 : Initialize CRSSANT for TG assembly

```python gapmcluster.py x.bedpe x_gapm.sam```

## OPTIONAL : Separate reads tagged by identifier

```awk '$0~/^@/ || $1~/-TAG/' x_pri_crssant.cliques.t_o0.2.sam > x_ crssant_TAG1.sam```

<p align="right">(<a href="#top">back to top</a>)</p>

# Ribosomal RNA (rRNA) DG asssembly

The abundance of rRNA reads yields useful data for comparative and quality control purposes. However, processing every rRNA associated read requires excessive computational resrouces. Thus, extracting a randomized subset of the reads for further analysis is a convenient workaround.

## STEP 1 : Extract rRNA reads

```
awk '$3=="hs45S"' x_prigap1_filtered.sam > x_prigap1_hs45S.sam
cut -f1-19 x_prigap1_hs45S.sam > x_prigap1_hs45S_tmp.sam
samtools view -H x_prigap1_filtered.sam > header
cat header x_prigap1_hs45S_tmp.sam > x_hs45S.sam
rm -f header x_prigap1_hs45S_tmp.sam
```

## STEP 2 : Extract a number of randomized reads

```python extract_random.py x_hs45S.sam @ 30000 x_hs45S_30000.sam```

## STEP 3 : Initialize CRSSANT for DG and NG assembly

```
samtools view -bS -o x_pri.bam x_hs45S_30000.sam
bedtools genomecov -bg -split -strand + -ibam x_hs45S_30000.bam -g staridxPath/chrNameLength.txt > x_plus.bedgraph
bedtools genomecov -bg -split -strand - -ibam x_hs45S_30000.bam -g staridxPath/chrNameLength.txt > x_minus.bedgraph
python crssant.py -cluster cliques -t_o 0.2 -out ./ x_hs45S_30000.sam Bed x_plus.bedgraph,x_minus.bedgraph
```

## STEP 4 : SAM to sorted BAM for IGV visualization

```
samtools view -bS -o x_pri.bam x_pri.sam
samtools sort x_pri.bam x_pri_sorted.bam
samtools index x_pri_sorted.bam
```

## STEP 5 : DG base pairs to stemloop arcs

```python dg_bp.py ref bedpe outname y/n DG_number```

## STEP 6 : Read coverage and DG arm endpoints

```python3 dg_endpoint.py input.sam```

## STEP 7 : Determine read arm nucleotide frequency

```python nt_freq.py sam ref + - rRNA all_reads reads ```

## STEP 8 : Determine read arm SHAPE reactivity frequency

```python SHAPE_freq.py inputsam shape_bedgraph DG_reads_cutoff DGcommon_ratio seqlen extendlen DG/reads outputprefix```

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

[Zhipeng Lu](https://zhipenglulab.org)

Project Link: [https://github.com/whl-usc/rna2d3d](https://github.com/whl-usc/rna2d3d)

<p align="right">(<a href="#top">back to top</a>)</p>
