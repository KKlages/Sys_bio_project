#!/bin/bash

# Loop over all subdirectories in the current directory
for dir in */; do
    echo "Checking $dir"
    
    # Check for .gz files one level below
    gz_files=$(find "$dir" -maxdepth 1 -type f -name "*.fastq.gz")
    
    if [ -n "$gz_files" ]; then
        echo "  Found .gz files in $dir"

        # Loop over all *_1.fastq.gz files
        for file in "$dir"/*_1.fastq.gz; do
            [ -e "$file" ] || continue  # Skip if no match

            # Extract prefix before the first underscore
            base=$(basename "$file")
            prefix=${base%%_*}
            
            # Create target directory
            target_dir="${dir}${prefix}"
            mkdir -p "$target_dir"
            
            # Move both _1 and _2 files
            for suffix in 1 2; do
                f="${dir}${prefix}_${suffix}.fastq.gz"
                if [ -f "$f" ]; then
                    mv "$f" "$target_dir/"
                    echo "    Moved $f → $target_dir/"
                fi
            done
        done
    fi
done

