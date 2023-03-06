#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=16G

module load gcc/11.3.0
module load openjdk 
module load perl/5.36.0

####################################################################################################
# contact:    wlee9829@gmail.com
# date:       2022_01_14
# script:     pre_process.sh
#    
# Shell script to standardize pre-processing PARIS- and SHARC-based sequencing data. Trimmomatic-0.36 
# and perl scripts from the icSHAPE pipeline are used. Change paths and variables accordingly.
#
# Trimmomatic: https://github.com/usadellab/Trimmomatic
# icSHAPE: https://github.com/qczhang/icSHAPE/tree/master/scripts
####################################################################################################

####################################################################################################
# Specify the following variables with absolute paths for software, working path, fasta/fastq.
####################################################################################################
# ScriptPath:	path to software (e.g., /bin/)
# WorkPath:	path for writing output files  (e.g., /data/preprocessing/test)
# fasta:	path to adapter.fa, contains adapter sequences (e.g., /preprocessing/P6SolexaRC35.fa)
# fastq:	path to "x".fastq.gz, the raw data file (e.g., /preprocessing/PARIS.fastq.gz)
# Outprefix:	prefix for output file names, "x" (e.g., PARIS)
####################################################################################################
ScriptPath=/bin
WorkPath=/data/preprocessing/
fasta=/ref/fasta/P6SolexaRC35.fa
fastq=/data/raw_data/L001_R1.fastq
Outprefix=test
####################################################################################################

####################################################################################################
#  Define the pre-processing function with preset parameters.
####################################################################################################

function preprocess() {
	mkdir $WorkPath
	cd $WorkPath

# Trimmomatic - remove 3'-end adapters
	java -jar $ScriptPath/trimmomatic/trimmomatic-0.36.jar SE -threads 16 -phred33 $fastq $Outprefix'_trim3.fastq' ILLUMINACLIP:$fasta:3:20:10 SLIDINGWINDOW:4:20 MINLEN:18
	echo

# splitFastq.pl - separate libraries based on barcode. PARIS and SHARC have barcodes R701-R724. Change these according to sequencing barcodes used. 
	echo "Splitting by barcodes."
	perl $ScriptPath/icSHAPE/icSHAPE-master/scripts/splitFastq.pl -U $Outprefix'_trim3.fastq' -b 6:6 -l \
	CGTGAT:R701::ACATCG:R702::GCCTAA:R703::TGGTCA:R704::CACTGT:R705::ATTGGC:R706::GATCTG:R707::TCAAGT:R708::CTGATC:R709::AAGCTA:R710::GTAGCC:R711::TACAAG:R712::TTGACT:R713::GGAACT:R714::TGACAT:R715::GGACGG:R716::CTCTAC:R717::GCGGAC:R718::TTTCAC:R719::GGCCAC:R720::CGAAAC:R721::CGTACG:R722::CCACTC:R723::GCTACC:R724 -s stats.txt -d ./
	echo

# readCollapse.pl - collapse the fastq file and remove any PCR duplicates
	for f in R*.fastq; do
		bcnn="${f%.fastq}"
		perl $ScriptPath/icSHAPE/icSHAPE-master/bin/readCollapse ${f} ${bcnn}'_trim3_nodup.fastq' 
		echo
	done

# Trimmomatic - remove 5'-head adapters
	for f in R*_trim3_nodup.fastq; do
	    bcnn="${f%.fastq}"
	    java -jar $ScriptPath/trimmomatic/trimmomatic-0.36.jar SE -threads 16 -phred33 ${f} $Outprefix'_trim_nodup_'${bcnn}'.fastq' HEADCROP:17 MINLEN:20
		echo
	done
}
####################################################################################################
preprocess
echo "Preprocessing completed."
####################################################################################################
