#!/bin/bash

echo "***** HPC job info ***** "
echo "Job ID: $LSB_JOBID"
echo "Job index within array: $LSB_JOBINDEX"
echo "Node: $(hostname)"
echo "Queue: $LSB_QUEUE"
echo "Job name: $LSB_JOBNAME"
echo "User: $LSB_USER"
echo "Submit directory: $LS_SUBCWD"
echo "Submission host: $LSB_SUB_HOST"
echo "Execution start time: $(date)"

# load modules
module load qtltools
module load tabix

# move to output directory
mkdir -p /rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/results/cis_eqtl

cd /rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/results/cis_eqtl

chr=$LSB_JOBINDEX
annot=$1
measure=$2
quant=$3
normalize=$4

prefix=$quant.${annot}.${measure}.chr${chr}
bed_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass2/files_for_analysis/bed_files/param_row_reads_1/${quant}_${annot}_gene_counts_${measure}.bed"
geno_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass1/files_for_analysis/1KG_vcf/genos_1kg_eur_500_snps_maf_0.01_chr${chr}.vcf.gz"

# bgzip and index the file if the gz file doesn't exist
if [[ ! -f "${bed_file}.gz" ]]; then
    echo "Compressing and indexing bed file"
    bgzip "$bed_file" && tabix -p bed "${bed_file}.gz"
else
    echo "tabix index already exists, skipping bgzip and tabix"
fi

if [[ "$normalize" == "T" ]]; then
    echo "Running normalization and QTLtools permutation at gene-level..."
    QTLtools cis \
        --vcf ${geno_file} \
        --bed ${bed_file}.gz \
        --permute 1000 \
        --out gene_qtls_perm_normY_${quant}_${prefix}.txt \
        --normal \
        --seed 1219
else
    QTLtools cis \
        --vcf ${geno_file} \
        --bed ${bed_file}.gz \
        --permute 1000 \
        --out gene_qtls_perm_normN_${quant}_${prefix}.txt \
        --seed 1219
fi


# QTLtools cis --vcf ${geno_file} --bed ${bed_file}.gz --permute 1000 --chunk 0 10 --out header_0.txt --normal --seed 1219

# unload modules
module unload tabix
module unload qtltools

echo "***** job ends ***** "
date

