#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G

module load gcc/11.3.0
module load openjdk
module load perl/5.36.0

####################################################################################################
# contact:    wlee9829@gmail.com
# date:       2023_03_06
# script:     staridx.sh
#
# This script standardizes the making of a star index for proximity-ligation and crosslinking data.
# STAR mapper is used, with function "genomeGenerate". Change paths and variables accordingly.
#
# STAR: https://github.com/alexdobin/STAR
####################################################################################################

####################################################################################################
# Specify the following variables with absolute paths for software, working path, fastq/gtf.
####################################################################################################
# ScriptPath:   path to software (e.g., /bin/)
# WorkPath:     path for writing output files  (e.g., /ref/staridx/)
# Fasta:        path to fasta file (e.g., /ref/fasta/)
####################################################################################################
ScriptPath=/bin/
WorkPath=/ref/staridx
Fasta=hg38mask14add.fa
####################################################################################################

####################################################################################################
#  Define the crssant function with preset parameters.
####################################################################################################
function index() {
        mkdir $WorkPath
        cd $WorkPath

        ## STAR genomeGenerate with parameters.
        $ScriptPath/STAR-2.7.1a/bin/Linux_x86_64/STAR \
        --runThreadN 8 \
        --runMode genomeGenerate \
        --genomeDir $WorkPath \
        --genomeFastaFiles $Fasta \
        --genomeSAindexNbases 14 \
        --genomeChrBinNbits 18 \
}
####################################################################################################
index
echo "StarIndex reference genome completed."
echo
####################################################################################################