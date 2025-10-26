#!/bin/bash

# RNA-seq processing pipeline: Download -> Convert -> QC -> Quantify

# Check required files and tools
if [ ! -f "samples.txt" ]; then
    echo "Error: samples.txt not found"
    exit 1
fi

if [ ! -f "crichardii.idx" ]; then
    echo "Error: crichardii.idx not found"
    exit 1
fi

# Check required tools
for tool in prefetch fastq-dump fastp kallisto; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool not found"
        exit 1
    fi
done

# Create output directory for final results
mkdir -p kallisto_results_final

# Process each sample
while IFS= read -r accession; do
    # Skip empty lines and comments
    if [[ -n "$accession" && ! "$accession" =~ ^[[:space:]]*# ]]; then
        echo "========================================="
        echo "Processing: $accession"
        echo "========================================="
        
        # Step 1: Download SRA file
        echo "Step 1: Downloading $accession..."
        prefetch "$accession"
        
        if [ ! -d "$accession" ]; then
            echo "Error: Failed to download $accession"
            continue
        fi
        
        # Navigate to sample directory
        cd "$accession"
        SRA_FILE="${accession}.sra"
        
        # Step 2: Convert SRA to FASTQ
        echo "Step 2: Converting SRA to FASTQ..."
        fastq-dump --split-files --gzip "$SRA_FILE"
        
        # Check if FASTQ files were created
        if [ ! -f "${accession}_1.fastq.gz" ] || [ ! -f "${accession}_2.fastq.gz" ]; then
            echo "Error: FASTQ conversion failed for $accession"
            cd ..
            rm -rf "$accession"
            continue
        fi
        
        # Step 3: Quality control and trimming
        echo "Step 3: Quality control and trimming..."
        fastp -i "${accession}_1.fastq.gz" -I "${accession}_2.fastq.gz" \
              -o "${accession}_1.clean.fastq.gz" -O "${accession}_2.clean.fastq.gz" \
              --length_required 50 \
              --html "fastp_report.html" \
              --json "fastp_report.json"
        
        # Check if cleaned files were created
        if [ ! -f "${accession}_1.clean.fastq.gz" ] || [ ! -f "${accession}_2.clean.fastq.gz" ]; then
            echo "Error: Quality trimming failed for $accession"
            cd ..
            rm -rf "$accession"
            continue
        fi
        
        # Step 4: Quantification with Kallisto
        echo "Step 4: Quantifying with Kallisto..."
        kallisto quant -i "../crichardii.idx" \
                      -o "kallisto_results" \
                      --bootstrap-samples=100 \
                      --threads=8 \
                      "${accession}_1.clean.fastq.gz" "${accession}_2.clean.fastq.gz"
        
        # Check if quantification succeeded
        if [ -d "kallisto_results" ] && [ -f "kallisto_results/abundance.tsv" ]; then
            echo "✓ Successfully processed $accession"
            
            # Copy results to final directory
            cp -r "kallisto_results" "../kallisto_results_final/${accession}_kallisto"
            
            # Go back to parent directory
            cd ..
            
            # Clean up: remove everything except kallisto results (which we've already copied)
            rm -rf "$accession"
            
        else
            echo "✗ Kallisto quantification failed for $accession"
            cd ..
            rm -rf "$accession"
        fi
        
        echo "Completed processing: $accession"
        echo ""
    fi
done < "samples.txt"

echo "========================================="
echo "Pipeline completed!"
echo "All results saved in: kallisto_results_final/"
echo "========================================="
