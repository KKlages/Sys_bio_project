#!/usr/bin/env bash
# SCRIPT: 03_extract_and_combine_hits.sh
# PURPOSE: Extract hit sequences from proteomes and combine them into a single file.

set -e -o pipefail

# --- Configuration ---
FASTA_DIR="data/01_protein_sequences/fern_proteomes"
RESULTS_DIR="results/phylogenetics/blast_hits"
OUTPUT_DIR="results/phylogenetics/hit_sequences"
COMBINED_FILE="results/phylogenetics/all_spl_candidates.fasta"

# --- Workflow ---
echo "Creating directory for individual hit sequences at ${OUTPUT_DIR}..."
mkdir -p "${OUTPUT_DIR}"

# 1. Run the Python script to extract sequences for each species
echo "Step 1: Extracting hit sequences using the Python helper script..."
python utils/select_hit_sequences.py --fasta_folder "${FASTA_DIR}" \
                                     --results_folder "${RESULTS_DIR}" \
                                     --output_folder "${OUTPUT_DIR}"

# 2. Combine all extracted hit sequences into one file (Missing Step)
echo "Step 2: Combining all hit sequences into a single file..."
# Clear the file if it already exists
> "${COMBINED_FILE}"
for file in "${OUTPUT_DIR}"/*.fa; do
  cat "$file" >> "${COMBINED_FILE}"
done

echo "--- Hit extraction and combination complete. ---"
echo "All candidate sequences are in: ${COMBINED_FILE}"