#!/bin/bash
set -euo pipefail

FOLDERS=("Ore" "Ala" "Aop")
for folder in "${FOLDERS[@]}"; do
    if [ -d "$folder" ]; then
        echo "Submitting SLURM job for $folder"
        sbatch submit_quant.slurm "$folder"
    else
        echo "Warning: Folder $folder not found, skipping."
    fi
done

