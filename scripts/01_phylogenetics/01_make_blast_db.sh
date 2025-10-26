#!/usr/bin/env bash
# SCRIPT: 01_make_blast_db.sh
# PURPOSE: Create a BLAST database from the Arabidopsis thaliana reference SPL proteins.

set -e -o pipefail

# --- Configuration ---
# Assumes your reference sequences are in the data directory
REF_SEQS="data/01_protein_sequences/athaliana_spl_known.fasta"
DB_DIR="results/phylogenetics/blast_db"
DB_NAME="${DB_DIR}/Ath_SPL_db"

# --- Workflow ---
echo "Creating BLAST database directory at ${DB_DIR}..."
mkdir -p "${DB_DIR}"

echo "Creating protein BLAST database from ${REF_SEQS}..."
makeblastdb -in "${REF_SEQS}" -dbtype prot -out "${DB_NAME}"

echo "--- BLAST database created successfully at ${DB_NAME} ---"