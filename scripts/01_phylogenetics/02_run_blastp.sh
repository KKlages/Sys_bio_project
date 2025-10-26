#!/usr/bin/env bash
# SCRIPT: 02_run_blastp.sh
# PURPOSE: Run BLASTp for each fern proteome against the A. thaliana SPL database.

set -e -o pipefail

# --- Configuration ---
# Assumes all target proteomes are in a single directory
INPUT_DIR="data/01_protein_sequences/fern_proteomes" 
DB_NAME="results/phylogenetics/blast_db/Ath_SPL_db"
OUTPUT_DIR="results/phylogenetics/blast_hits"
E_VALUE="1e-7"

# --- Workflow ---
echo "Creating output directory for BLAST results at ${OUTPUT_DIR}..."
mkdir -p "${OUTPUT_DIR}"

echo "Running BLASTp for all FASTA files in ${INPUT_DIR}..."

for file in "${INPUT_DIR}"/*.fa; do
  # Extract the base name of the file (e.g., "Azolla_filiculoides")
  base_name=$(basename "${file}" .fa)
  output_file="${OUTPUT_DIR}/${base_name}_results.tsv"
  
  echo "  -> Processing ${base_name}..."
  
  blastp -query "$file" \
         -db "${DB_NAME}" \
         -out "${output_file}" \
         -evalue "${E_VALUE}" \
         -outfmt 6 \
         -max_target_seqs 1 # As per your original script
done

echo "--- BLASTp search complete. Results are in ${OUTPUT_DIR} ---"