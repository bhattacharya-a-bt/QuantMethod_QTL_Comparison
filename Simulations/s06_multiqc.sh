#!/bin/bash

module load multiqc
module load R/4.3.1

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate multiqc-1.13

PASS=$1
PARAM_ROW_READS=$2

out_dir=/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/fastqc/param_row_reads_${PARAM_ROW_READS}
cd $out_dir

multiqc .

Rscript /rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/scripts/s06_multiqc.R ${PASS} $PARAM_ROW_READS

module unload multiqc
module unload R

conda deactivate