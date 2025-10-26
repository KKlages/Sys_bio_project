# SCRIPT: helper_filter_cdsearch_results.py
# PURPOSE: Merges and filters hitdata files from NCBI CD-Search.
#
# It performs a stringent filtering step:
# 1. Merges all 'hitdata-part*.txt' files in the directory.
# 2. Identifies any protein query that has AT LEAST ONE domain hit
#    that is NOT 'SBP' or 'SBP superfamily'.
# 3. Removes ALL entries for those identified proteins.
# 4. Outputs a final, clean list of protein IDs that only ever hit SBP domains.
#
# REQUIREMENTS: pandas (pip install pandas)
#
# USAGE:
# python helper_filter_cdsearch_results.py -o <output_id_list.txt>
#
# EXAMPLE:
# python helper_filter_cdsearch_results.py -o ../../results/phylogenetics/final_spl_protein_ids.txt

import pandas as pd
import glob
import argparse
import os

def filter_cdsearch_results(output_file: str):
    """
    Loads all 'hitdata-part*.txt' files, filters them, and saves a list of clean query IDs.
    
    Args:
        output_file (str): Path to save the final list of clean protein IDs.
    """
    # Find all hitdata files in the current directory
    hitdata_files = glob.glob("hitdata*.txt")
    if not hitdata_files:
        print("Error: No 'hitdata*.txt' files found in this directory.")
        print("Please download hit details from NCBI CD-Search and place them here.")
        return

    print(f"Found files: {', '.join(hitdata_files)}")
    
    all_dfs = []
    for file_path in hitdata_files:
        try:
            # Dynamically find the header row, as NCBI can add comment lines
            with open(file_path, 'r', errors='ignore') as f:
                for i, line in enumerate(f):
                    if line.startswith("Query\tHit type"):
                        header_line = i
                        break
                else: # no-break
                    print(f"Warning: Header row 'Query\tHit type' not found in {file_path}. Skipping.")
                    continue
            
            df = pd.read_csv(file_path, sep='\t', skiprows=header_line)
            all_dfs.append(df)
        except Exception as e:
            print(f"Could not process {file_path}: {e}")
    
    if not all_dfs:
        print("Error: No data could be loaded from the hitdata files. Exiting.")
        return

    merged_df = pd.concat(all_dfs, ignore_index=True)

    # --- Filtering Logic ---
    # 1. Define the only domains we consider valid for an SPL protein
    allowed_domains = ["SBP", "SBP superfamily"]

    # 2. Find all 'Query' values that are associated with a non-allowed domain
    invalid_queries = merged_df.loc[
        ~merged_df["Short name"].isin(allowed_domains), "Query"
    ].unique()

    # 3. Filter the main dataframe to remove ALL rows belonging to these invalid queries
    cleaned_df = merged_df[~merged_df["Query"].isin(invalid_queries)]
    
    # 4. Get the final list of unique, clean query IDs.
    # The ID is the first part of the 'Query' string, e.g., ">protein_id ..."
    unique_clean_ids = cleaned_df["Query"].str.split(n=1).str[0].str.lstrip('>').unique()
    
    # 5. Save the list to the specified output file
    try:
        with open(output_file, "w") as f:
            for query_id in sorted(unique_clean_ids):
                f.write(query_id + "\n")
                
        print("\n--- Filtering Summary ---")
        print(f"Total unique queries in raw data: {merged_df['Query'].nunique()}")
        print(f"Queries removed due to non-SBP domains: {len(invalid_queries)}")
        print(f"Final number of clean SPL queries: {len(unique_clean_ids)}")
        print(f"Clean ID list saved to: {os.path.abspath(output_file)}")

    except Exception as e:
        print(f"Error writing to output file {output_file}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter NCBI CD-Search results to get high-confidence SPL protein IDs.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-o", "--output_file",
        required=True,
        help="Path to save the final list of clean protein IDs (e.g., final_ids.txt)."
    )
    
    args = parser.parse_args()
    filter_cdsearch_results(args.output_file)