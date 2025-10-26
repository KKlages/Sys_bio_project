
#!/bin/bash

# Loop over all batch CSV files
for batch in batch_*.csv; do
    echo "[$(date)] Starting downloads for $batch"
    ./download_batch.sh "$batch"

    # Optional: Check if all expected FASTQ files exist (basic check)
    missing=0
    awk -F',' 'NR>1{print $2}' "$batch" | while read -r cmd; do
        # Extract output filename from wget command
        fname=$(echo "$cmd" | awk '{print $NF}' | sed "s/'//g")
        if [[ ! -f "$fname" ]]; then
            echo "WARNING: Expected file $fname not found after download for $batch"
            missing=1
        fi
    done
    if [[ $missing -eq 1 ]]; then
        echo "ERROR: Some files missing after download for $batch. Skipping RSEM for this batch."
        continue
    fi

    echo "[$(date)] Downloads finished for $batch. Submitting RSEM job..."
    sbatch rsem_batch.slurm "$batch"
    echo "[$(date)] RSEM job submitted for $batch"
done

echo "All batches processed."

