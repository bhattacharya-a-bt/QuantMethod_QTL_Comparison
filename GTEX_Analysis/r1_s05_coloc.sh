#!/bin/bash

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

# assign analysis parameters

PSPACE="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/requant_paramspace.txt"

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"
echo "pheno_file=$1"
echo "pheno_name=$2"
echo "liftoverTo38=$3"
echo "n_bins=$4"
echo "bin_index=$5"

echo "Running colocalization script..."
Rscript /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/Scripts/r1_s05_coloc.R ${annot} ${quant} ${tissue} $1 $2 $3 $4 $5

# unload modules
module unload R
module unload qtltools
module unload tabix

echo "***** job ends ***** "
date

