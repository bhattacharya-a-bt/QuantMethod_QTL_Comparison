#!/bin/bash

# Usage: ./count_gtex_bams.sh <tissue> <user>

# ============================ USER CONFIGURATION ============================
# Edit the defaults below, or export these variables before running.
#   sample_IDs  GTEx v8 sample attributes file (tab-delimited; column 1 is the
#               sample ID, column 3 is the tissue)
#   bamdir      directory holding the source GTEx BAM files
# ============================================================================
sample_IDs="${sample_IDs:-/path/to/GTEx_v8/GTEx_v8_sample_attributes.txt}"
bamdir="${bamdir:-/path/to/GTEx/SourceFiles/Bam}"

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
