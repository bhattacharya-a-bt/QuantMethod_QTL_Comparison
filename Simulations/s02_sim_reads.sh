#!/bin/bash

# ============================ USER CONFIGURATION ============================
# Edit the default below, or export SCRIPT_DIR before submitting the job.
#   SCRIPT_DIR  directory containing s02_sim_reads.R
# ============================================================================
SCRIPT_DIR="${SCRIPT_DIR:-/path/to/QuantMethod_QTL_Comparison/Simulations}"

# load modules
module load R/4.3.1

PASS=$1
PARAM_ROW_READS=$2
BASE_DIR=$3
TXOME_FA=$4

# run the simulationR script with the job index as an argument
Rscript ${SCRIPT_DIR}/s02_sim_reads.R ${LSB_JOBINDEX} ${PASS} ${PARAM_ROW_READS} ${BASE_DIR} ${TXOME_FA}

# unload modules
module unload R

echo "***** job ends ***** "
date
