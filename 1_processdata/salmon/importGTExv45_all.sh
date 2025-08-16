#!/bin/bash
#BSUB -J gtex_reprocess[1-48]
#BSUB -o /rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess/logs/bam_to_fastq_%J_%I.out
#BSUB -e /rsrch5/scratch/epi/abhattacharya3/GTEx_reprocess/logs/bam_to_fastq_%J_%I.err
#BSUB -W 2:00
#BSUB -n 10
#BSUB -M 80
#BSUB -R rusage[mem=80]
#BSUB -q short

# Load required modules (if needed)
module load R

# Run the R script with the current array index
Rscript /rsrch5/home/epi/abhattacharya3/processGTEx/importGTExv45.R $LSB_JOBINDEX
