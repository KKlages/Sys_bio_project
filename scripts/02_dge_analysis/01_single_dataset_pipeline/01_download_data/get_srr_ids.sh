#!/bin/bash

# List of SRR IDs
srr_ids=(
SRR6480271
SRR6480272
SRR6480273
SRR6480274
SRR6480275
SRR6480276
SRR6480277
SRR6480278
SRR6480279
SRR6480280
SRR6480281
SRR6480282
)

# Download each SRR using prefetch
for srr in "${srr_ids[@]}"; do
    echo "Downloading $srr..."
    prefetch "$srr"
done

echo "All downloads complete."

