#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=epyc-64
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --job-name=GEO2fastq

###############################################################################
###############################################################################
## Use this script to quickly download data from a GEO study link.            #
##                                                                            #
## From any GEO link associated with dataset, open GEO Accession viewer.      #
## Select SRA Run Selector, download the Total: Accession List                #
##                                                                            #
## SRR_LIST_FILE: The text file containing SRR IDs                            #
## OUTPUTDIR: The output PATH for the downloaded data                         #
###############################################################################
## Change the specified paths below accordingly:
SRR_LIST_FILE="SRR_Acc_List.txt"
OUTPUT_DIR="./geo_rnaseq_data"
###############################################################################
###############################################################################

# Function to check if data is single-end or paired-end using fastq-dump
check_sequencing_type() {
    local SRR_ID=$1
    local SRA_PATH="${OUTPUT_DIR}/${SRR_ID}/${SRR_ID}.sra"
    
    # Check if the SRA file exists
    if [[ ! -f "${SRA_PATH}" ]]; then
        echo "SRA file not found: ${SRA_PATH}"
        return 1
    fi
    
    # Extract a small subset of reads using fastq-dump
    local read_data=$(fastq-dump -X 1 --split-spot --stdout "${SRA_PATH}")
    
    # Determine sequencing type from the output
    if echo "${read_data}" | grep -q "@.*/1"; then
        echo "PAIRED"
    elif echo "${read_data}" | grep -q "@"; then
        echo "SINGLE"
    else
        echo "UNKNOWN"
    fi
}

# Check if SRA Toolkit is installed
if ! command -v prefetch &> /dev/null || ! command -v fasterq-dump &> /dev/null || ! command -v vdb-dump &> /dev/null
then
    echo "SRA Toolkit could not be found. Please install it first."
    exit 1
fi

# Check if SRR_Acce_List.txt exists
if [[ ! -f "$SRR_LIST_FILE" ]]; then
    echo "File $SRR_LIST_FILE not found!"
    exit 1
fi

# Create output directory if it doesn't already exist
mkdir -p "$OUTPUT_DIR"

# Read each SRR ID from the file and process 
while IFS= read -r SRR_ID; do

    if [[ ! -z "$SRR_ID" ]]; then
        echo "Downloading $SRR_ID"
        
        if ! prefetch "$SRR_ID" -O "$OUTPUT_DIR"; then
            echo "Failed to download $SRR_ID. Skipping..."
            continue
        fi

        echo "Converting $SRR_ID to FASTQ"

        # Determine if the data is single-end or paired-end
        sequencing_type=$(check_sequencing_type "$SRR_ID")
        if [[ "$sequencing_type" == "PAIRED" ]]; then
            # Paired-end data, split into separate files
            if ! fastq-dump --split-files -O "$OUTPUT_DIR" --gzip "${OUTPUT_DIR}/${SRR_ID}/${SRR_ID}.sra"; then
                echo "Failed to split $SRR_ID to separate files. Skipping..."
                continue
            fi
        elif [[ "$sequencing_type" == "SINGLE" ]]; then
            # Single-end data, no need to split files
            if ! fastq-dump -O "$OUTPUT_DIR" --gzip "${OUTPUT_DIR}/${SRR_ID}/${SRR_ID}.sra"; then
                echo "Failed to convert $SRR_ID to FASTQ. Skipping..."
                continue
            fi
        else
            echo "Unable to determine sequencing type for $SRR_ID. Skipping..."
            continue
        fi

        # Optionally remove the .sra file to save space
        rm -r "${OUTPUT_DIR}/${SRR_ID}"
    fi
done < "$SRR_LIST_FILE"

echo "All downloads and conversions are completed at: $(date +%H:%M:%S)."