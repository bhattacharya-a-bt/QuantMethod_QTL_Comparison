#!/bin/bash

# Usage: ./cleanup_GR.sh <tissue> <user>
#
# WARNING: this permanently deletes the per-sample intermediates for a tissue.
# Only run it once that tissue's aggregated RDS files exist.

# ============================ USER CONFIGURATION ============================
# Edit the defaults below, or export these variables before running.
#   SCRATCH_DIR  per-user scratch working directory to clear
#   REQUANT_DIR  requantification output root; the raw per-sample quant
#                subdirectories under each annotation are removed
# ============================================================================
tissue="$1"
user="$2"

SCRATCH_DIR="${SCRATCH_DIR:-/path/to/scratch/${user}/GTEx_v8}"
REQUANT_DIR="${REQUANT_DIR:-/path/to/GTEx_v8/requants}"

rm -rf ${SCRATCH_DIR}/temp/${tissue}
rm -rf ${SCRATCH_DIR}/raw/${tissue}
rm -rf ${SCRATCH_DIR}/reports/${tissue}
rm -rf ${SCRATCH_DIR}/metrics/${tissue}
rm -rf ${SCRATCH_DIR}/align/${tissue}
rm -rf ${SCRATCH_DIR}/quant/${tissue}
rm -rf ${REQUANT_DIR}/GENCODE_v27/${tissue}
rm -rf ${REQUANT_DIR}/GENCODE_v38/${tissue}
rm -rf ${REQUANT_DIR}/GENCODE_v45/${tissue}
rm -rf ${REQUANT_DIR}/Ensembl/${tissue}