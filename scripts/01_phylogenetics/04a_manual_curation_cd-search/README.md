# Manual Curation using NCBI Conserved Domain Search (CD-Search)

This directory contains helpers for the manual but critical step of validating SPL candidates. The goal is to remove any protein that has a predicted domain other than the SBP domain.

### Workflow

**Step 1: Split Candidate File for Web Submission**

The NCBI CD-Search web tool has a limit on the number of sequences you can submit at once. Use the `helper_split_fasta.py` script to break down your large candidate file (`results/phylogenetics/all_spl_candidates.fasta`) into smaller chunks.

```bash
# Example: Split the file into 3 parts
python helper_split_fasta.py \
  --input_file ../../results/phylogenetics/all_spl_candidates.fasta \
  --num_files 3 \
  --output_prefix candidates_part
This will create candidates_part_1.fa, candidates_part_2.fa, etc.
Step 2: Submit to NCBI CD-Search
Go to the NCBI CD-Search website.
For each candidates_part_*.fa file, upload it and run the search.
After each search completes, click the "Download" tab and download the "Hit Details (text)" file.
Rename these files descriptively (e.g., hitdata-part1.txt, hitdata-part2.txt) and save them in this directory.
Step 3: Filter the CD-Search Results
Once you have all your hitdata-part*.txt files, use the helper_filter_cdsearch_results.py script. It will merge all hit data files, identify any protein that has a non-SBP domain associated with it, and generate a final list of clean protein IDs.
code
Bash
# The script will automatically find all hitdata-part*.txt files in the directory
python helper_filter_cdsearch_results.py --output_file ../../results/phylogenetics/final_spl_protein_ids.txt
This final .txt file contains the list of high-confidence SPL protein IDs.
Step 4: Create Final FASTA File
Use the generated final_spl_protein_ids.txt to extract the corresponding sequences from your all_spl_candidates.fasta file. You can use a simple grep tool like seqkit.
code
Bash
# Run from the 01_phylogenetics/ directory
seqkit grep -f ../results/phylogenetics/final_spl_protein_ids.txt \
           ../results/phylogenetics/all_spl_candidates.fasta \
           -o ../results/phylogenetics/final_curated_spls_unaligned.fasta

# IMPORTANT: Don't forget to add your original Arabidopsis sequences back in for the alignment!
cat ../data/01_protein_sequences/athaliana_spl_known.fasta >> ../results/phylogenetics/final_curated_spls_unaligned.fasta