#!/bin/bash
#BSUB -J "clean_qtl_[1-22]"
#BSUB -o /path/to/scratch/GTEx_gencode_comp/pass2/logs/qtl/get_true_beta_info_%I.out
#BSUB -e /path/to/scratch/GTEx_gencode_comp/pass2/logs/qtl/get_true_beta_info_%I.err
#BSUB -n 1
#BSUB -M 30
#BSUB -R "rusage[mem=30G]"
#BSUB -W 24:00
#BSUB -q medium

# ============================ USER CONFIGURATION ============================
# Edit the default below, or export SCRIPT_DIR before submitting the job.
#   SCRIPT_DIR  directory containing s10_clean_qtl_res.R
# The #BSUB -o/-e log paths above must be edited directly: LSF does not expand
# shell variables in #BSUB directives.
# s10_clean_qtl_res.R has its own configuration block for its paths.
# ============================================================================
SCRIPT_DIR="${SCRIPT_DIR:-/path/to/QuantMethod_QTL_Comparison/Simulations}"

# load modules
module load R/4.3.1

CHR=$LSB_JOBINDEX

# run the simulationR script with the job index as an argument
Rscript ${SCRIPT_DIR}/s10_clean_qtl_res.R ${CHR}

# unload modules
module unload R

echo "***** job ends ***** "
date
