#!/bin/bash
set -e

# Define base and analysis directory
BASE_DIR="/mnt/ceph-hdd/projects/scc_ubmg_devries/KKlages/mrna_subproject"
ANALYSIS_DIR="$BASE_DIR/analysis"

# Create main analysis directory
mkdir -p "$ANALYSIS_DIR"

echo "🔍 Searching species folders..."

# Loop through species directories (e.g., Aob, Adi, etc.)
for species_dir in "$BASE_DIR"/*/; do
    species_name=$(basename "$species_dir")

    # Skip non-directory entries or unrelated files
    [[ ! -d "$species_dir/kallisto_quant" ]] && continue

    echo "Processing $species_name..."

    # Create species-specific subfolder inside analysis/
    target_dir="$ANALYSIS_DIR/$species_name"
    mkdir -p "$target_dir"

    # Find all .tsv files under kallisto_quant/ERR*/ and copy them
    find "$species_dir/kallisto_quant" -type f -name "*.tsv" -exec cp {} "$target_dir/" \;

    echo "Copied .tsv files for $species_name."
done

echo "All .tsv files collected under $ANALYSIS_DIR/"

