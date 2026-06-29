#!/bin/bash
#BSUB -J "qtl_perm[337-352]"
#BSUB -o /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/logs/qtl/qtl_perm_%J_%I.out
#BSUB -e /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/logs/qtl/qtl_perm_%J_%I.out
#BSUB -q medium
#BSUB -W 24:00
#BSUB -n 1
#BSUB -M 30
#BSUB -R rusage[mem=30]

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
module load R/4.3.1
module load tabix

# assign analysis parameters

# need to make sure run tabix on bed file

PSPACE="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/requant_paramspace.txt"

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"

# move to output directory
mkdir -p /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/cis_eqtl_results

cd /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/cis_eqtl_results

prefix=${tissue}.${annot}.${quant}
bed_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/${tissue}/${annot}_${quant}.v8.normalized_expression.bed"
cov_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/requant_analyses/${tissue}/${tissue}_formatted_covariates.txt"
geno_file="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/GTEx_838_v8_maf0.01_autosomes_unrelated.vcf.gz" # copy from GTEx_MAF_0.01_passQC

# bgzip and index the file if the gz file doesn't exist
if [[ ! -f "${bed_file}.gz" ]]; then
    echo "Compressing and indexing bed file"
    bgzip "$bed_file" && tabix -p bed "${bed_file}.gz"
else
    echo "tabix index already exists, skipping bgzip and tabix"
fi

echo "Running normalization and QTLtools permutation at gene-level..."
QTLtools cis --vcf ${geno_file} --bed ${bed_file}.gz --permute 1000 --cov ${cov_file} --out gene_qtls_perm_normalized_${prefix}.txt --normal --seed 1219

# unload modules
module unload R
module unload tabix
module unload qtltools

echo "***** job ends ***** "
date

