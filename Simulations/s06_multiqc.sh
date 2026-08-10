#!/bin/bash

# ============================ USER CONFIGURATION ============================
# Edit the defaults below, or export these variables before submitting the job.
#   BASE_DIR    base project directory; the pass number is appended internally
#   SCRIPT_DIR  directory containing s06_multiqc.R
#   CONDA_BIN   path to the conda binary
# s06_multiqc.R has its own configuration block for its paths.
# ============================================================================
BASE_DIR="${BASE_DIR:-/path/to/scratch/GTEx_gencode_comp}"
SCRIPT_DIR="${SCRIPT_DIR:-/path/to/QuantMethod_QTL_Comparison/Simulations}"
CONDA_BIN="${CONDA_BIN:-/path/to/miniforge3/bin/conda}"

module load multiqc
module load R/4.3.1

eval "$(${CONDA_BIN} shell.bash hook)"
conda activate multiqc-1.13

PASS=$1
PARAM_ROW_READS=$2

out_dir=${BASE_DIR}/pass${PASS}/files_for_analysis/fastqc/param_row_reads_${PARAM_ROW_READS}
cd $out_dir

multiqc .

Rscript ${SCRIPT_DIR}/s06_multiqc.R ${PASS} $PARAM_ROW_READS

module unload multiqc
module unload R

conda deactivate