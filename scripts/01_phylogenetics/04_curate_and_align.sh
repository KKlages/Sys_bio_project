#!/usr/bin/env bash
# SCRIPT: 04_curate_and_align.sh
# PURPOSE: Document data curation and run multiple sequence alignment with MAFFT.

set -e -o pipefail

# --- Configuration ---
CANDIDATE_SEQS="results/phylogenetics/all_spl_candidates.fasta"
# THIS IS THE KEY MANUAL STEP: You need to create this file after curation.
FINAL_SEQS="results/phylogenetics/final_curated_spls.fasta"
ALIGNED_SEQS="results/phylogenetics/aligned_spl_proteins.fa"

# --- Workflow ---
echo "--- Starting Curation and Alignment ---"
echo "IMPORTANT MANUAL STEP:"
echo "The file '${CANDIDATE_SEQS}' contains all raw BLAST hits."
echo "You must now manually curate this file:"
echo "  1. Add the original A. thaliana reference sequences."
echo "  2. Use your 'remove_duplicate_ids.py' script or other methods to clean the data."
echo "  3. Perform NCBI CD-search to validate the SBP domain."
echo "  4. Save the final, clean set of sequences to '${FINAL_SEQS}'."

# Check if the user has created the final file
if [ ! -f "$FINAL_SEQS" ]; then
    echo -e "\nError: Final sequence file not found at '${FINAL_SEQS}'. Please perform the manual curation steps above and then re-run this script."
    exit 1
fi

echo -e "\nFinal sequence file found. Proceeding with alignment..."

# Run MAFFT for multiple sequence alignment
echo "Running MAFFT for alignment..."
mafft --auto --thread -1 "${FINAL_SEQS}" > "${ALIGNED_SEQS}"

echo "--- Alignment complete. Output file: ${ALIGNED_SEQS} ---"