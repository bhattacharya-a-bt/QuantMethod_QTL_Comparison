#!/bin/bash

# Usage: ./loop_process_short_reads_GR.sh <tissue> <user>

# ============================ USER CONFIGURATION ============================
# Edit the default below, or export DIR_SCRIPTS before running.
#   DIR_SCRIPTS  directory containing autoLoop_process_short_reads_GR.lsf
# ============================================================================
DIR_SCRIPTS="${DIR_SCRIPTS:-/path/to/QuantMethod_QTL_Comparison/GTEX_Analysis/Requantification_Scripts}"

tissue="$1"
user="$2"

cd ${DIR_SCRIPTS}

bsub -env tissue=${tissue},user=${user} < autoLoop_process_short_reads_GR.lsf
