# SCRIPT: helper_split_fasta.py
# PURPOSE: Splits a large multi-FASTA file into a specified number of smaller files.
#
# This is useful for submitting sequences to web servers (like NCBI CD-Search)
# that have a limit on input file size or sequence count.
#
# REQUIREMENTS: BioPython (pip install biopython)
#
# USAGE:
# python helper_split_fasta.py -i <input_file.fa> -n <num_files> -o <output_prefix>
#
# EXAMPLE:
# python helper_split_fasta.py -i ../../results/phylogenetics/all_spl_candidates.fasta -n 3 -o candidates_part

import argparse
import math
from Bio import SeqIO

def split_fasta(input_file: str, num_files: int, output_prefix: str):
    """
    Splits a FASTA file into a specified number of smaller files.

    Args:
        input_file (str): Path to the input FASTA file.
        num_files (int): The number of output files to create.
        output_prefix (str): The prefix for the output filenames.
                             (e.g., 'part' -> 'part_1.fa', 'part_2.fa')
    """
    try:
        # First pass: count the number of sequences to determine split size
        print(f"Counting sequences in {input_file}...")
        all_records = list(SeqIO.parse(input_file, "fasta"))
        total_records = len(all_records)
        if total_records == 0:
            print("Input file contains no sequences. Exiting.")
            return

        seqs_per_file = math.ceil(total_records / num_files)
        print(f"Found {total_records} total sequences.")
        print(f"Splitting into {num_files} files with ~{seqs_per_file} sequences each.")

        # Second pass: write the records to the new files
        file_count = 1
        record_index = 0
        for i in range(0, total_records, seqs_per_file):
            chunk = all_records[i:i + seqs_per_file]
            output_filename = f"{output_prefix}_{file_count}.fa"
            
            with open(output_filename, "w") as out_handle:
                SeqIO.write(chunk, out_handle, "fasta")
            
            print(f"Wrote {len(chunk)} sequences to {output_filename}")
            file_count += 1
            
        print("\nSplitting complete.")

    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_file}'")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a multi-FASTA file into smaller chunks.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i", "--input_file",
        required=True,
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "-n", "--num_files",
        required=True,
        type=int,
        help="The number of smaller files to create."
    )
    parser.add_argument(
        "-o", "--output_prefix",
        required=True,
        help="Prefix for the output files (e.g., 'output_part')."
    )

    args = parser.parse_args()

    if args.num_files <= 0:
        print("Error: Number of files must be a positive integer.")
    else:
        split_fasta(args.input_file, args.num_files, args.output_prefix)