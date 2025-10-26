#!/bin/bash
set -e

# Define base analysis folder (where renamed files are)
BASE_DIR="/mnt/ceph-hdd/projects/scc_ubmg_devries/KKlages/mrna_subproject/analysis"
OUTPUT_FILE="combined_metadata.tsv"

# Header
echo -e "species\tsample\ttissue\treplicate\tgroup" > "$OUTPUT_FILE"

# Python to walk through files and extract metadata
python3 <<EOF >> "$OUTPUT_FILE"
import os
import re

base_dir = "$BASE_DIR"

for species in sorted(os.listdir(base_dir)):
    species_dir = os.path.join(base_dir, species)
    if not os.path.isdir(species_dir):
        continue

    for filename in sorted(os.listdir(species_dir)):
        if not filename.endswith(".tsv"):
            continue

        match = re.match(r"(.+)_([0-9]+)\.tsv$", filename)
        if not match:
            print(f"Skipping unrecognized file: {filename}", file=os.sys.stderr)
            continue

        tissue, replicate = match.groups()
        sample = f"{tissue}_{replicate}"
        group = tissue  # or change this if you want custom groups
        print(f"{species}\t{sample}\t{tissue}\t{replicate}\t{group}")
EOF

echo "Metadata saved to $OUTPUT_FILE"

