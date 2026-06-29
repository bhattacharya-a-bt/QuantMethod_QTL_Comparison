#!/bin/bash

# Usage: ./count_gtex_bams.sh <tissue> <user>

sample_IDs=/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx_v8_sample_attributes.txt
bamdir=/rsrch3/scratch/reflib/GTEx/SourceFiles/Bam

tissue="$1"

# Collect sample IDs for this tissue
FILES=( $(awk -F'\t' -v t="$tissue" 'NR>1 && $3 == t {print $1}' "${sample_IDs}") )

count=0

for FILE in "${FILES[@]}"; do
    # Find a matching BAM (same pattern your original pipeline used)
    MATCH=$(ls "${bamdir}/${FILE}"*.bam 2>/dev/null | head -n 1)
    
    if [[ -n "$MATCH" ]]; then
        ((count++))
    fi
done

echo "Tissue: ${tissue}"
echo "Samples with matching BAMs: ${count}"
echo "Total expected sample IDs: ${#FILES[@]}"
