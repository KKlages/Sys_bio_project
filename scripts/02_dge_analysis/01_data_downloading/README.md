# Data Downloading Scripts

This directory contains various scripts used to download raw sequencing data.

-   **`get_ids.txt`**: A bash script that uses `prefetch` to download SRA files based on a hardcoded list of SRR accessions.
-   **`split.py`**: A Python script that reads a master CSV file of download links (`download_links.csv`) and splits it into smaller `batch_*.csv` files, one for each species abbreviation.
-   **`download_batch.txt`**: A shell script designed to be called by `loop.txt`. It takes a batch CSV file as input and runs the `wget` commands listed in it in parallel using `xargs`.
-   **`loop.txt`**: The main script to orchestrate the batch download process. It loops through all `batch_*.csv` files and calls `download_batch.txt` for each one.