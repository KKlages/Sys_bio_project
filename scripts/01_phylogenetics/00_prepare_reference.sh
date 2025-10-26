#!/usr/bin/env bash
# SCRIPT: 00_prepare_reference.sh
# PURPOSE: Extract A. thaliana SPL protein sequences from a full proteome FASTA file.

set -e -o pipefail

# --- Configuration ---
# Assumes the full Araport11 proteome is in the data directory
FULL_PROTEOME="data/01_protein_sequences/Araport11_pep_20250411.fasta"
OUTPUT_FASTA="data/01_protein_sequences/athaliana_spl_known.fasta"

# --- Prerequisite Check ---
if [ ! -f "$FULL_PROTEOME" ]; then
    echo "Error: Full proteome file not found at ${FULL_PROTEOME}"
    exit 1
fi

echo "Extracting representative SPL gene models from ${FULL_PROTEOME}..."

# --- Workflow ---
# This awk command is taken from your notebook. It's the most reliable method.
# It finds headers containing "Symbols:" and "SPL[number]" and prints the entire FASTA record.
awk 'BEGIN { RS=">"; ORS="" } /Symbols:.*SPL[0-9]+/ { print ">" $0 }' "${FULL_PROTEOME}" > "${OUTPUT_FASTA}"

# Count the sequences to confirm
NUM_SEQS=$(grep -c ">" "${OUTPUT_FASTA}")
echo "--- Reference preparation complete. ---"
echo "Found ${NUM_SEQS} SPL sequences."
echo "Reference file saved to: ${OUTPUT_FASTA}"