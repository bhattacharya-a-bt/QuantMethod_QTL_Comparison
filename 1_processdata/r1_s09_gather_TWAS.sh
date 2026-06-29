#!/bin/bash
#BSUB -J "gatherTWAS[1-768]"
#BSUB -o /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/logs/twas/gatherTWAS_%J_%I.out
#BSUB -e /rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/logs/twas/gatherTWAS_%J_%I.err
#BSUB -q short
#BSUB -W 3:00
#BSUB -n 1
#BSUB -M 10
#BSUB -R rusage[mem=10]

echo "***** HPC job info ***** "
echo "Job ID: $LSB_JOBID"
echo "Job index within array: $LSB_JOBINDEX"
echo "Node: $(hostname)"
echo "Queue: $LSB_QUEUE"
echo "Job name: $LSB_JOBNAME"
echo "Submit directory: $LS_SUBCWD"
echo "Submission host: $LSB_SUB_HOST"
echo "Execution start time: $(date)"

# assign analysis parameters

PSPACE="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/requant_paramspace.txt"
pheno="Height"

read -r annot quant tissue < <(awk -v row="$LSB_JOBINDEX" 'NR==row {print $1, $2, $3}' ${PSPACE})

echo "annot=$annot"
echo "quant=$quant"
echo "tissue=$tissue"
echo "pheno=$pheno"

INDIR="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/twas_results/${tissue}/${annot}"
OUTDIR="/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/requants/twas_results"
OUTFILE="${OUTDIR}/r1_twas_z_${annot}_${quant}_${tissue}_${pheno}.txt"

echo "Concatenating files..."

# concatenate all matching files in this subfolder
find "$INDIR" \
  -type f \
  -name "*_${pheno}_${annot}_${quant}.txt.gz" \
  -print0 | sort -z | \
xargs -0 zcat > "$OUTFILE"

echo "***** job ends ***** "
date

