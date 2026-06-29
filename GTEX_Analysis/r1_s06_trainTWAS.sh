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
module load qtltools
module unload R
module load R/4.3.1

# assign analysis parameters

PSPACE="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/requant_paramspace.txt"

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"

echo "Running TWAS model training..."
Rscript /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/Scripts/r1_s06_trainTWAS.R ${annot} ${quant} ${tissue} $1

# unload modules
module unload R
module unload qtltools

echo "***** job ends ***** "
date

