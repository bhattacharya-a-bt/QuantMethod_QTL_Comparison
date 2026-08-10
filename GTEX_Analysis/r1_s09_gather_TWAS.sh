#!/bin/bash
#BSUB -J "gatherTWAS[1-768]"
#BSUB -o /path/to/GTEx_v8/requants/logs/twas/gatherTWAS_%J_%I.out
#BSUB -e /path/to/GTEx_v8/requants/logs/twas/gatherTWAS_%J_%I.err
#BSUB -q short
#BSUB -W 3:00
#BSUB -n 1
#BSUB -M 10
#BSUB -R rusage[mem=10]

echo "***** HPC job info ***** "
echo "Job ID: $LSB_JOBID"
echo "Job index within array: $LSB_JOBINDEX"
echo "Node: $(hostname)"
echo "Queue: $LSB_QUEUE"
echo "Job name: $LSB_JOBNAME"
echo "Submit directory: $LS_SUBCWD"
echo "Submission host: $LSB_SUB_HOST"
echo "Execution start time: $(date)"

# ============================ USER CONFIGURATION ============================
# Edit the defaults below, or export these variables before submitting the job.
#   PSPACE        parameter space file (columns: annot quant tissue), one row per
#                 array index
#   TWAS_OUT_DIR  TWAS results directory written by r1_s08_runTWAS.R; per-gene
#                 files live in '<tissue>/<annot>/' and the concatenated file is
#                 written at the top level
#   pheno         phenotype label to gather
# The #BSUB -o/-e log paths above must be edited directly: LSF does not expand
# shell variables in #BSUB directives.
# ============================================================================
PSPACE="${PSPACE:-/path/to/GTEx_v8/requants/requant_paramspace.txt}"
TWAS_OUT_DIR="${TWAS_OUT_DIR:-/path/to/GTEx_v8/requants/twas_results}"
pheno="${pheno:-Height}"

# assign analysis parameters

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"
echo "pheno=$pheno"

INDIR="${TWAS_OUT_DIR}/${tissue}/${annot}"
OUTDIR="${TWAS_OUT_DIR}"
OUTFILE="${OUTDIR}/r1_twas_z_${annot}_${quant}_${tissue}_${pheno}.txt"

echo "Concatenating files..."

# concatenate all matching files in this subfolder
find "$INDIR" \
  -type f \
  -name "*_${pheno}_${annot}_${quant}.txt.gz" \
  -print0 | sort -z | \
xargs -0 zcat > "$OUTFILE"

echo "***** job ends ***** "
date

