#!/bin/bash
#BSUB -J "make_bed[273-288]"
#BSUB -o /path/to/GTEx_v8/requants/logs/qtl/make_bed_%J_%I.out
#BSUB -e /path/to/GTEx_v8/requants/logs/qtl/make_bed_%J_%I.out
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1
#BSUB -M 20
#BSUB -R rusage[mem=20]

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
#   PSPACE      parameter space file (columns: annot quant tissue), one row per
#               array index
#   SCRIPT_DIR  directory containing r1_s02_RDStoBed.R
# The #BSUB -o/-e log paths above must be edited directly: LSF does not expand
# shell variables in #BSUB directives.
# r1_s02_RDStoBed.R has its own configuration block for its input/output paths.
# ============================================================================
PSPACE="${PSPACE:-/path/to/GTEx_v8/requants/requant_paramspace.txt}"
SCRIPT_DIR="${SCRIPT_DIR:-/path/to/QuantMethod_QTL_Comparison/GTEX_Analysis}"

# load modules
module load R/4.3.1

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"

Rscript ${SCRIPT_DIR}/r1_s02_RDStoBed.R ${annot} ${quant} ${tissue}

# unload modules
module unload R

echo "***** job ends ***** "
date
