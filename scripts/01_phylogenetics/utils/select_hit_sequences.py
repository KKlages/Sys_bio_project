# SCRIPT: select_hit_sequences.py
# PURPOSE: Extracts FASTA sequences that are listed as hits in BLAST result files.
#
# This script iterates through all BLAST result files (e.g., ending in '_results.tsv')
# in a specified directory. For each result file, it finds the corresponding original
# FASTA proteome file, reads the list of hit IDs, and writes the full FASTA records
# for those hits to a new file in an output directory.
#
# REQUIREMENTS: BioPython (pip install biopython)
#
# USAGE:
# python select_hit_sequences.py \
#   --fasta_folder <path_to_proteomes> \
#   --results_folder <path_to_blast_results> \
#   --output_folder <path_to_save_hit_sequences>
#
# EXAMPLE (run from scripts/01_phylogenetics/ directory):
# python utils/select_hit_sequences.py \
#   --fasta_folder ../../data/01_protein_sequences/fern_proteomes \
#   --results_folder ../../results/phylogenetics/blast_hits \
#   --output_folder ../../results/phylogenetics/hit_sequences

import os
import glob
import argparse
from Bio import SeqIO

def extract_hit_sequences(fasta_folder: str, results_folder: str, output_folder: str):
    """
    Extracts sequences from FASTA files based on BLAST results.

    Args:
        fasta_folder (str): Path to the directory containing the original proteome FASTA files.
        results_folder (str): Path to the directory containing BLAST result files (outfmt 6).
        output_folder (str): Path to the directory where extracted sequences will be saved.
    """
    # Create the output folder if it doesn't exist
    try:
        os.makedirs(output_folder, exist_ok=True)
    except OSError as e:
        print(f"Error creating output directory {output_folder}: {e}")
        return

    # Use a consistent file extension for blast results, e.g., .tsv
    blast_results_pattern = os.path.join(results_folder, "*_results.tsv")
    result_files = glob.glob(blast_results_pattern)

    if not result_files:
        print(f"Warning: No BLAST result files found matching '{blast_results_pattern}'")
        return

    print(f"Found {len(result_files)} BLAST result files to process.")

    # Process each results file
    for results_file in result_files:
        try:
            # Extract the base name (e.g., "Azolla_filiculoides") from the results file name
            base_name = os.path.basename(results_file).replace("_results.tsv", "")
            
            # Construct the path to the corresponding original FASTA file
            fasta_file = os.path.join(fasta_folder, f"{base_name}.fa")
            
            # Skip if the corresponding FASTA file doesn't exist
            if not os.path.exists(fasta_file):
                print(f"Warning: No matching FASTA file found for {results_file}. Looked for '{fasta_file}'. Skipping.")
                continue
            
            # Read the target IDs from the results file (second column of outfmt 6)
            hit_ids = set()
            with open(results_file, 'r') as f:
                for line in f:
                    if line.strip():  # Skip empty lines
                        # BLAST outfmt 6 format is: qseqid sseqid ...
                        # The hit in the database is the second column (sseqid)
                        hit_ids.add(line.strip().split('\t')[1])
            
            if not hit_ids:
                print(f"No hits found in {results_file}. Skipping.")
                continue

            # Create the path for the output file
            output_file = os.path.join(output_folder, f"{base_name}_hit_sequences.fa")
            
            # Find matching sequences and write them to the output file
            hits_found_in_fasta = 0
            with open(output_file, 'w') as out_f:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    # The record.id is the part of the header before the first space
                    if record.id in hit_ids:
                        SeqIO.write(record, out_f, "fasta")
                        hits_found_in_fasta += 1
            
            print(f"Processed {base_name}: Found {hits_found_in_fasta}/{len(hit_ids)} requested IDs and saved to {os.path.basename(output_file)}")

        except Exception as e:
            print(f"An error occurred while processing {os.path.basename(results_file)}: {e}")

    print("\nExtraction process complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences based on BLAST hit IDs.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--fasta_folder",
        required=True,
        help="Directory containing the original proteome FASTA files (e.g., species.fa)."
    )
    parser.add_argument(
        "--results_folder",
        required=True,
        help="Directory containing the BLAST output files (e.g., species_results.tsv)."
    )
    parser.add_argument(
        "--output_folder",
        required=True,
        help="Directory where the output FASTA files with hit sequences will be saved."
    )

    args = parser.parse_args()
    
    extract_hit_sequences(args.fasta_folder, args.results_folder, args.output_folder)