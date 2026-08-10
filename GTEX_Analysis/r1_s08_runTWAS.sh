#!/bin/bash

# ============================ USER CONFIGURATION ============================
# Edit the default below, or export SCRIPT_DIR before submitting the job.
#   SCRIPT_DIR  directory containing r1_s08_runTWAS.R
# r1_s08_runTWAS.R has its own configuration block for its input/output paths.
# ============================================================================
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
module load tabix
module load qtltools
module unload R
module load R/4.3.1

# analysis parameters

echo "index=$LSB_JOBINDEX"
echo "pheno_file=$1"
echo "pheno_name=$2"
echo "liftoverTo38=$3"
echo "n_bins=$4"
echo "bin_index=$5"

echo "Running TWAS..."
Rscript ${SCRIPT_DIR}/r1_s08_runTWAS.R $LSB_JOBINDEX $1 $2 $3 $4 $5

# unload modules
module unload R
module unload qtltools
module unload tabix

echo "***** job ends ***** "
date

