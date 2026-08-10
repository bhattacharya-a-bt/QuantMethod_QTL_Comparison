#!/bin/bash

# Usage: ./run_quants_GR.sh <tissue> <user>

# ============================ USER CONFIGURATION ============================
# Edit the default below, or export DIR_SCRIPTS before running.
#   DIR_SCRIPTS  directory containing run_quants_GR.lsf
# ============================================================================
DIR_SCRIPTS="${DIR_SCRIPTS:-/path/to/QuantMethod_QTL_Comparison/GTEX_Analysis/Requantification_Scripts}"

tissue="$1"
user="$2"

bsub -env TISSUE=${tissue},USER=${user} < ${DIR_SCRIPTS}/run_quants_GR.lsf
