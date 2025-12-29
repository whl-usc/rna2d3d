#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=epyc-64
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=GEO2fastq

###############################################################################
###############################################################################
## Use this script to download data from a GEO study via SRA accessions.      ##
##                                                                           ##
## From the GEO Accession Viewer:                                            ##
##   → SRA Run Selector                                                      ##
##   → Download “Total: Accession List”                                      ##
##                                                                           ##
## SRR_LIST_FILE : Text file containing SRR IDs (one per line)               ##
## OUTPUT_DIR   : Output directory for FASTQ files                           ##
###############################################################################
###############################################################################

# ----------------------------- USER SETTINGS --------------------------------

SRR_LIST_FILE="SRR_Acc_List.txt"
OUTPUT_DIR="./geo_rnaseq_data"
THREADS=${SLURM_CPUS_PER_TASK:-8}

# -----------------------------------------------------------------------------

set -o pipefail

# ----------------------------- SANITY CHECKS ---------------------------------

for cmd in prefetch fasterq-dump gzip; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "ERROR: $cmd not found. Please load/install SRA Toolkit."
        exit 1
    fi
done

if [[ ! -f "$SRR_LIST_FILE" ]]; then
    echo "ERROR: SRR list file not found: $SRR_LIST_FILE"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ----------------------------- MAIN LOOP -------------------------------------

while IFS= read -r SRR_ID; do
    [[ -z "$SRR_ID" ]] && continue

    echo "============================================================"
    echo "Processing $SRR_ID"
    echo "Start time: $(date)"

    # --------------------------- DOWNLOAD ------------------------------------

    echo "Downloading SRA file for $SRR_ID"

    if ! prefetch --max-size 0 "$SRR_ID" -O "$OUTPUT_DIR"; then
        echo "ERROR: Failed to download $SRR_ID"
        continue
    fi

    SRA_PATH="${OUTPUT_DIR}/${SRR_ID}/${SRR_ID}.sra"

    if [[ ! -f "$SRA_PATH" ]]; then
        echo "ERROR: SRA file not found after download: $SRA_PATH"
        continue
    fi

    # --------------------------- CONVERSION ----------------------------------

    echo "Converting $SRR_ID to FASTQ"

    if ! fasterq-dump \
        --split-files \
        --threads "$THREADS" \
        --outdir "$OUTPUT_DIR" \
        "$SRA_PATH"; then
        echo "ERROR: FASTQ conversion failed for $SRR_ID"
        continue
    fi

    # --------------------------- COMPRESSION ---------------------------------

    echo "Compressing FASTQ files for $SRR_ID"
    gzip "${OUTPUT_DIR}/${SRR_ID}"*.fastq

    # --------------------------- CLEANUP -------------------------------------

    echo "Removing SRA files for $SRR_ID"
    rm -rf "${OUTPUT_DIR:?}/${SRR_ID}"

    echo "Finished $SRR_ID at $(date)"
    echo "============================================================"

done < "$SRR_LIST_FILE"

echo "All downloads and conversions completed at $(date)"