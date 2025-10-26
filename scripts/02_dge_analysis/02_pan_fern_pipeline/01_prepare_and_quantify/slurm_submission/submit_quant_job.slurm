#!/bin/bash
#SBATCH --job-name=rsem_batch
#SBATCH --output=rsem_batch_%j.out
#SBATCH --error=rsem_batch_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=medium


source activate fern_rnaseq_env

FOLDER=$1
echo "Starting job in folder: $FOLDER"
cd "$FOLDER" || { echo "ERROR: Cannot enter folder $FOLDER"; exit 1; }

bash ./run_quant.sh

echo "Job finished in folder: $FOLDER"

