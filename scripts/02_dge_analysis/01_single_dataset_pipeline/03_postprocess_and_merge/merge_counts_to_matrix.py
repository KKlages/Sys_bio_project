import pandas as pd
import glob
import os

# Set path to your abundance files
path = ""
files = sorted(glob.glob(os.path.join(path, "*.tsv")))

# Initialize an empty DataFrame
all_counts = pd.DataFrame()

# Loop over each file
for f in files:
    sample_name = os.path.basename(f).replace(".tsv", "")
    df = pd.read_csv(f, sep="\t", usecols=["target_id", "est_counts"])
    df.rename(columns={"est_counts": sample_name}, inplace=True)

    if all_counts.empty:
        all_counts = df
    else:
        all_counts = all_counts.merge(df, on="target_id")

# Set target_id as index
all_counts.set_index("target_id", inplace=True)

# Save to file for edgeR
output_file = os.path.join(path, "counts.txt")
all_counts.to_csv(output_file, sep="\t")