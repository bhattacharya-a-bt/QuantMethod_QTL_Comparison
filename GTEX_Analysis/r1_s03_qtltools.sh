#!/bin/bash
#BSUB -J "qtl_perm[337-352]"
#BSUB -o /path/to/GTEx_v8/requants/logs/qtl/qtl_perm_%J_%I.out
#BSUB -e /path/to/GTEx_v8/requants/logs/qtl/qtl_perm_%J_%I.out
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

# ============================ USER CONFIGURATION ============================
# Edit the defaults below, or export these variables before submitting the job.
#   PSPACE        parameter space file (columns: annot quant tissue), one row per
#                 array index
#   ANALYSIS_DIR  directory holding the per-tissue BED and covariate files
#                 written by r1_s02_RDStoBed.R and r1_s01_writeCovariates.R
#   EQTL_OUT_DIR  output directory for the QTLtools permutation results
#   GENO_FILE     bgzipped+tabixed GTEx genotype VCF
#                 (copy from GTEx_MAF_0.01_passQC)
# The #BSUB -o/-e log paths above must be edited directly: LSF does not expand
# shell variables in #BSUB directives.
# ============================================================================
PSPACE="${PSPACE:-/path/to/GTEx_v8/requants/requant_paramspace.txt}"
ANALYSIS_DIR="${ANALYSIS_DIR:-/path/to/scratch/GTEx_gencode_comp/requant_analyses}"
EQTL_OUT_DIR="${EQTL_OUT_DIR:-/path/to/GTEx_v8/requants/cis_eqtl_results}"
GENO_FILE="${GENO_FILE:-/path/to/scratch/GTEx_gencode_comp/GTEx_838_v8_maf0.01_autosomes_unrelated.vcf.gz}"

# load modules
module load qtltools
module load R/4.3.1
module load tabix

# assign analysis parameters

# need to make sure run tabix on bed file

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"

# move to output directory
mkdir -p ${EQTL_OUT_DIR}

cd ${EQTL_OUT_DIR}

prefix=${tissue}.${annot}.${quant}
bed_file="${ANALYSIS_DIR}/${tissue}/${annot}_${quant}.v8.normalized_expression.bed"
cov_file="${ANALYSIS_DIR}/${tissue}/${tissue}_formatted_covariates.txt"
geno_file="${GENO_FILE}"

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

