#!/bin/bash
BATCH_CSV="$1"
LOG="download_$(basename "$BATCH_CSV" .csv).log"

# Extract wget commands from column 2 (skip header) and run 2 in parallel
awk -F',' 'NR>1{print $2}' "$BATCH_CSV" | \
    xargs -I CMD -P 2 bash -c "echo [\$(date)] CMD | tee -a '$LOG'; CMD >>'$LOG' 2>&1"

