#!/bin/bash

# Exit script on any error
set -e
# Treat unset variables as an error
set -u
# Ensure pipeline errors are propagated
set -o pipefail

# --- Configuration ---
# Working directory (where SRR directories are located)
PROJECT_DIR="$(pwd)"
KALLISTO_INDEX_FILENAME="azolla.idx"
THREADS=8
BASE_RESULTS_DIR_NAME="${PROJECT_DIR}/results"
# --- End Configuration ---

# Construct full paths
KALLISTO_INDEX_PATH="${PROJECT_DIR}/${KALLISTO_INDEX_FILENAME}"
BASE_RESULTS_DIR="${BASE_RESULTS_DIR_NAME}"

# Get the directory where the script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Define SRR samples (directories that should already exist)
SRR_SAMPLES=(
    "SRR6480271"
    "SRR6480272"
    "SRR6480273"
    "SRR6480274"
    "SRR6480275"
    "SRR6480276"
    "SRR6480277"
    "SRR6480278"
    "SRR6480279"
    "SRR6480280"
    "SRR6480281"
    "SRR6480282"
)

# --- Sanity Checks ---
if [ ! -f "$KALLISTO_INDEX_PATH" ]; then
    echo "ERROR: Kallisto index not found at ${KALLISTO_INDEX_PATH}"
    echo "Please ensure '${KALLISTO_INDEX_FILENAME}' exists at the specified path."
    exit 1
fi

# Check for required tools
for cmd in fastq-dump fastp kallisto; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "ERROR: Required command '$cmd' not found in PATH."
        exit 1
    fi
done

echo "RNA-Seq Processing Pipeline Started"
echo "-----------------------------------"
echo "Script directory: ${SCRIPT_DIR}"
echo "Project directory: ${PROJECT_DIR}"
echo "Kallisto index: ${KALLISTO_INDEX_PATH}"
echo "Base results directory: ${BASE_RESULTS_DIR}"
echo "Threads: ${THREADS}"
echo "-----------------------------------"

# Create the base results directory if it doesn't exist
mkdir -p "${BASE_RESULTS_DIR}"

# Process each SRR sample
for SRR_ID in "${SRR_SAMPLES[@]}"; do
    
    echo ""
    echo "=================================================="
    echo "Processing SRR ID: ${SRR_ID}"
    echo "=================================================="

    # Check if SRR directory exists
    SRR_INPUT_DIR="${PROJECT_DIR}/${SRR_ID}"
    if [ ! -d "${SRR_INPUT_DIR}" ]; then
        echo "ERROR: SRR directory not found: ${SRR_INPUT_DIR}. Skipping this sample."
        continue
    fi

    # Check if SRA file exists
    SRA_FILE="${SRR_INPUT_DIR}/${SRR_ID}.sra"
    if [ ! -f "${SRA_FILE}" ]; then
        echo "ERROR: SRA file not found: ${SRA_FILE}. Skipping this sample."
        continue
    fi

    # --- Define paths for the current SRR ID ---
    SRR_SPECIFIC_DIR="${BASE_RESULTS_DIR}/${SRR_ID}"
    FASTQ_RAW_DIR="${SRR_SPECIFIC_DIR}/fastq_raw"
    FASTQ_CLEAN_DIR="${SRR_SPECIFIC_DIR}/fastq_clean"
    KALLISTO_OUTPUT_DIR="${SRR_SPECIFIC_DIR}/kallisto_quant"

    # FASTQ files after SRA conversion (now gzipped)
    RAW_FASTQ_1="${FASTQ_RAW_DIR}/${SRR_ID}_1.fastq.gz"
    RAW_FASTQ_2="${FASTQ_RAW_DIR}/${SRR_ID}_2.fastq.gz"

    # Output FASTQ files from fastp
    CLEAN_FASTQ_1="${FASTQ_CLEAN_DIR}/${SRR_ID}_1.clean.fastq"
    CLEAN_FASTQ_2="${FASTQ_CLEAN_DIR}/${SRR_ID}_2.clean.fastq"
    FASTP_HTML_REPORT="${FASTQ_CLEAN_DIR}/${SRR_ID}_fastp.html"
    FASTP_JSON_REPORT="${FASTQ_CLEAN_DIR}/${SRR_ID}_fastp.json"

    # --- Create directories for the current SRR_ID ---
    echo "Creating directories for ${SRR_ID}..."
    mkdir -p "${FASTQ_RAW_DIR}"
    mkdir -p "${FASTQ_CLEAN_DIR}"
    mkdir -p "${KALLISTO_OUTPUT_DIR}"

    # 1. Convert SRA to FASTQ with gzip compression
    echo "-----------------------------------"
    echo "Step 1: Converting SRA to gzipped FASTQ for ${SRR_ID}..."
    if fastq-dump --split-files --gzip --outdir "${FASTQ_RAW_DIR}" "${SRA_FILE}"; then
        echo "SRA to gzipped FASTQ conversion complete for ${SRR_ID}"
    else
        echo "ERROR: SRA to FASTQ conversion failed for ${SRR_ID}. Skipping this sample."
        rm -rf "${SRR_SPECIFIC_DIR}"
        continue
    fi

    # Verify FASTQ files exist
    if [ ! -f "${RAW_FASTQ_1}" ] || [ ! -f "${RAW_FASTQ_2}" ]; then
        echo "ERROR: Converted FASTQ files not found for ${SRR_ID}. Skipping this sample."
        rm -rf "${SRR_SPECIFIC_DIR}"
        continue
    fi

    echo "Converted files:"
    echo "  R1: ${RAW_FASTQ_1}"
    echo "  R2: ${RAW_FASTQ_2}"

    # 2. Quality control and trimming
    echo "-----------------------------------"
    echo "Step 2: Quality control and trimming with fastp for ${SRR_ID}..."
    if fastp \
        -i "${RAW_FASTQ_1}" \
        -I "${RAW_FASTQ_2}" \
        -o "${CLEAN_FASTQ_1}" \
        -O "${CLEAN_FASTQ_2}" \
        --length_required 50 \
        --thread "${THREADS}" \
        --html "${FASTP_HTML_REPORT}" \
        --json "${FASTP_JSON_REPORT}"; then
        echo "fastp complete. Cleaned FASTQ in ${FASTQ_CLEAN_DIR}"
    else
        echo "ERROR: fastp failed for ${SRR_ID}. Skipping this sample."
        rm -rf "${SRR_SPECIFIC_DIR}"
        continue
    fi

    # 3. Transcript quantification
    echo "-----------------------------------"
    echo "Step 3: Quantifying transcripts with Kallisto for ${SRR_ID}..."
    if kallisto quant \
        -i "${KALLISTO_INDEX_PATH}" \
        -o "${KALLISTO_OUTPUT_DIR}" \
        --threads="${THREADS}" \
        --bootstrap-samples=100 \
        "${CLEAN_FASTQ_1}" \
        "${CLEAN_FASTQ_2}"; then
        echo "Kallisto quantification complete. Output in ${KALLISTO_OUTPUT_DIR}"
    else
        echo "ERROR: Kallisto quantification failed for ${SRR_ID}. Skipping this sample."
        rm -rf "${FASTQ_RAW_DIR}" "${FASTQ_CLEAN_DIR}"
        find "${KALLISTO_OUTPUT_DIR}" -empty -type d -delete
        find "${SRR_SPECIFIC_DIR}" -empty -type d -delete
        continue
    fi

    # 4. Clean up intermediate files
    echo "-----------------------------------"
    echo "Step 4: Cleaning up intermediate files for ${SRR_ID}..."
    rm -rf "${FASTQ_RAW_DIR}" "${FASTQ_CLEAN_DIR}"
    echo "Intermediate files removed. Kept: ${KALLISTO_OUTPUT_DIR}"

    echo "--------------------------------------------------"
    echo "Successfully processed SRR ID: ${SRR_ID}"
    echo "--------------------------------------------------"

done

echo ""
echo "==================================="
echo "All samples processed."
echo "Final Kallisto outputs are in subdirectories under: ${BASE_RESULTS_DIR}"
echo "==================================="