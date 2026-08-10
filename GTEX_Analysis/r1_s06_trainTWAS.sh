#!/bin/bash

# ============================ USER CONFIGURATION ============================
# Edit the defaults below, or export these variables before submitting the job.
#   PSPACE      parameter space file (columns: annot quant tissue), one row per
#               array index
#   SCRIPT_DIR  directory containing r1_s06_trainTWAS.R
# r1_s06_trainTWAS.R has its own configuration block for its paths.
# ============================================================================
PSPACE="${PSPACE:-/path/to/GTEx_v8/requants/requant_paramspace.txt}"
SCRIPT_DIR="${SCRIPT_DIR:-/path/to/QuantMethod_QTL_Comparison/GTEX_Analysis}"

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
module unload R
module load R/4.3.1

# assign analysis parameters

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"

echo "Running TWAS model training..."
Rscript ${SCRIPT_DIR}/r1_s06_trainTWAS.R ${annot} ${quant} ${tissue} $1

# unload modules
module unload R
module unload qtltools

echo "***** job ends ***** "
date

