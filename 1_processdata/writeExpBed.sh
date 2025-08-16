#!/bin/bash
#BSUB -J "writeExpBed[1-48]"
#BSUB -o /rsrch5/scratch/epi/abhattacharya3/GTEx_compare/logs/writeExpBed_%J_%I.out
#BSUB -e /rsrch5/scratch/epi/abhattacharya3/GTEx_compare/logs/writeExpBed_%J_%I.err
#BSUB -q medium
#BSUB -W 4:00
#BSUB -n 1
#BSUB -M 10
#BSUB -R rusage[mem=80]


# Load R module (adjust module name as needed for your system)
module load R plink tabix htslib

# Run the R script with the job array index
Rscript writeExpBed_v27.R ${LSB_JOBINDEX}

echo "Job ${LSB_JOBINDEX} completed successfully"