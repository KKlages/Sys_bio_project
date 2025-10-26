#!/bin/bash
set -e

# Base project path
BASE_DIR="/mnt/ceph-hdd/projects/scc_ubmg_devries/KKlages/mrna_subproject"
METADATA="sample_table.tsv"  # <- your TSV file here

# Python script to do the renaming across all species
python3 <<EOF
import os
import csv
from collections import defaultdict

base_dir = "$BASE_DIR"
metadata_file = "$METADATA"

# Normalize tissue names for filesystem
def normalize(tissue):
    return tissue.strip().replace(" ", "_").replace("/", "_")

# Read metadata with tab delimiter
tissue_counts = defaultdict(int)
err_to_info = {}

with open(metadata_file, newline='') as tsvfile:
    reader = csv.DictReader(tsvfile, delimiter='\t')
    for row in reader:
        err_id = row["ERR"]
        species_abbrev = row["Abbrev"]
        tissue = normalize(row["Tissue"])
        err_to_info[err_id] = (species_abbrev, tissue)

# Process each ERR ID
for err_id, (species_abbrev, tissue) in err_to_info.items():
    species_path = os.path.join(base_dir, species_abbrev)
    subdir = os.path.join(species_path, "kallisto_quant", err_id)
    src = os.path.join(subdir, "abundance.tsv")

    if not os.path.isfile(src):
        print(f"⚠️  Skipping missing file for {err_id}: {src}")
        continue

    tissue_counts[(species_abbrev, tissue)] += 1
    rep = tissue_counts[(species_abbrev, tissue)]
    new_name = f"{tissue}_{rep}.tsv"
    dst = os.path.join(subdir, new_name)

    print(f"✅ Renaming {src} -> {dst}")
    os.rename(src, dst)
EOF

echo "✅ All files renamed successfully."

