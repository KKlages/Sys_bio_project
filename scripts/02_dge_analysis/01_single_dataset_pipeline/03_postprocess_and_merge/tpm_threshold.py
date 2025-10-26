import pandas as pd

def threshold_tpm(input_csv, output_csv):
    # Read the CSV file
    df = pd.read_csv(input_csv)
    
    # Assume first column is 'id' or similar, do not threshold it
    id_col = df.columns[0]
    
    # Apply threshold to all other columns (set values < 1 to 0)
    df.iloc[:, 1:] = df.iloc[:, 1:].applymap(lambda x: 0 if x < 1 else x)
    
    # Save the result to a new CSV file
    df.to_csv(output_csv, index=False)
    print(f"Processed {input_csv} -> {output_csv}")

# Example usage for all three files:
threshold_tpm("Afi_PRJEB45223.csv", "Afi_PRJEB45223_thresh.csv")
threshold_tpm("Afi_PRJNA430527.csv", "Afi_PRJNA430527_thresh.csv")
threshold_tpm("Cri.csv", "Cri_thresh.csv")

