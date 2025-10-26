#!/bin/bash
#SBATCH --job-name=iqtree_spl_analysis
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/iqtree_%j.out
#SBATCH --error=logs/iqtree_%j.err

# SCRIPT: 05_run_iqtree.sh
# PURPOSE: Run IQ-TREE2 for phylogenetic analysis on an HPC using SLURM.
# USAGE: sbatch 05_run_iqtree.sh

# --- Setup ---
# Create a log directory if it doesn't exist
mkdir -p logs

# Activate your conda environment
# Replace 'fern_rnaseq_env' with the name from your environment.yml
source activate fern_rnaseq_env

# --- Configuration ---
ALIGNED_SEQS="results/phylogenetics/aligned_spl_proteins.fa"
MODEL="JTT+G4" # Determined by ModelFinder in your report
THREADS=16

# --- Analysis Options ---
# Choose which analysis to run by uncommenting one of the blocks below.

# OPTION 1: Standard model test and tree search (as in run_iqtree.txt)
# This is a good first run to confirm the model and get a basic tree.
# echo "Running standard IQ-TREE analysis..."
# iqtree2 -s "${ALIGNED_SEQS}" -m "${MODEL}" -T ${THREADS} --prefix results/phylogenetics/iqtree_standard

# OPTION 2: Ultrafast Bootstrap Analysis (as in ultrafast_tree.txt)
# The recommended analysis for getting robust branch supports.
# The -bnni flag optimizes branches after UFBoot, which is good practice.
echo "Running IQ-TREE with Ultrafast Bootstrap..."
iqtree2 -s "${ALIGNED_SEQS}" -m "${MODEL}" -T ${THREADS} -B 1000 -bnni --prefix results/phylogenetics/iqtree_ultrafast

# OPTION 3: Resume a failed job (as in un_iqtree_resume.txt)
# To use this, you must have a checkpoint file from a previous failed run.
# The --prefix must match the original run you want to resume.
# echo "Attempting to resume previous IQ-TREE run..."
# iqtree2 --resume --prefix results/phylogenetics/iqtree_ultrafast


echo "--- IQ-TREE job submitted. ---"