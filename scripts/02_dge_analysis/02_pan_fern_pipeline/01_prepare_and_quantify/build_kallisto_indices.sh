#!/bin/bash

for dir in */; do
  fasta=$(find "$dir" -maxdepth 1 -name "*.fasta")
  if [[ -n "$fasta" ]]; then
    base=$(basename "$fasta" .fasta)
    index_output="${dir%/}/${base}_kallisto.idx"
    echo "Indexing $fasta -> $index_output"
    kallisto index -i "$index_output" "$fasta"
  fi
done

