#!/bin/bash

# ============================ USER CONFIGURATION ============================
# Edit the default below, or export SCRIPT_DIR before submitting the job.
#   SCRIPT_DIR  directory containing s05_featureCounts.R
# s05_featureCounts.R has its own configuration block for its paths.
# ============================================================================
SCRIPT_DIR="${SCRIPT_DIR:-/path/to/QuantMethod_QTL_Comparison/Simulations}"

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
module load R/4.3.1

PASS=$1
GENO_PASS=$2
PARAM_ROW_READS=$3
GTF=$4
ANNOT_NAME=$5

Rscript ${SCRIPT_DIR}/s05_featureCounts.R ${PASS} ${GENO_PASS} $PARAM_ROW_READS $GTF $ANNOT_NAME

# unload modules
module unload R

echo "***** job ends ***** "
date