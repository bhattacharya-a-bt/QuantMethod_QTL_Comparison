#!/bin/bash

# ============================ USER CONFIGURATION ============================
# Edit the default below, or export SCRIPT_DIR before submitting the job.
#   SCRIPT_DIR  directory containing s06_export_tximeta.R
# s06_export_tximeta.R has its own configuration block for its paths.
# ============================================================================
SCRIPT_DIR="${SCRIPT_DIR:-/path/to/QuantMethod_QTL_Comparison/Simulations}"

# load modules
module load R/4.3.1

PASS=$1
PARAM_ROW_READS=$2
ANNOT=$3

# run the simulationR script with the job index as an argument
Rscript ${SCRIPT_DIR}/s06_export_tximeta.R ${PASS} ${PARAM_ROW_READS} ${ANNOT}

# unload modules
module unload R

echo "***** job ends ***** "
date
