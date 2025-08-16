# GTEx Colocalization Analysis Pipeline

This repository contains scripts for performing colocalization analysis between expression quantitative trait loci (eQTLs) from GTEx data and genome-wide association study (GWAS) signals for various cancer types.

## Overview

The pipeline identifies genes where eQTL signals colocalize with cancer GWAS signals, suggesting shared causal variants that affect both gene expression and cancer risk. The analysis uses eCAVIAR for colocalization testing across multiple GTEx versions and processing methods.

## Files Description

### Shell Scripts

#### `runQTLTools.sh`
- **Purpose**: LSF job submission script for running eQTL analysis using QTLtools
- **Job Array**: 1-48 jobs (one per tissue)
- **Resources**: 1 CPU, 80GB memory, 200-hour time limit
- **Queue**: long
- **Functions**:
  - Runs both permutation (1000 permutations) and nominal pass analysis
  - Processes all BED files (different GTEx versions/methods) within each tissue
  - Uses random seeds for reproducible permutation testing
- **Output**: `*_cisQTL.txt` (permutation results) and `*_cisQTL_nominal.txt` (nominal p-values)

#### `runColoc.sh`
- **Purpose**: LSF job submission script for running colocalization analysis
- **Job Array**: 1-1536 jobs (processes genes in batches of 25)
- **Resources**: 1 CPU, 40GB memory, 48-hour time limit
- **Queue**: long
- **Calls**: `runColoc.R` with array index

#### `countColoc.sh` 
- **Purpose**: LSF job submission script for summarizing colocalization results
- **Job Array**: 1-480 jobs (48 tissues × 10 cancers)
- **Resources**: 10 CPUs, 80GB memory, 24-hour time limit
- **Queue**: medium
- **Calls**: `countColoc.R` with tissue and cancer indices

### R Scripts

#### `runColoc.R`
**Main colocalization analysis script**

**Key Functions:**
- Processes genes in batches (25 genes per job)
- Performs genomic coordinate liftover (hg19 → hg38) for GWAS data
- Extracts linkage disequilibrium (LD) from GTEx reference panel
- Runs eCAVIAR colocalization analysis for each gene-GWAS-tissue combination

**Input Data:**
- **eGenes**: Pre-identified expression genes from `eGene_colocalization.RDS`
- **GWAS Files**: 10 cancer GWAS summary statistics:
  - Endometrial, Renal, Prostate, Ovarian, Melanoma
  - Lung, Head & Neck, Glioma, Colorectal, Breast
- **GTEx Data**: eQTL nominal p-values across 6 conditions:
  - 3 versions: v27, v38, v45
  - 2 methods: STAR, Salmon RNA-seq processing

**Analysis Steps:**
1. Load gene coordinates and GWAS data
2. Liftover GWAS coordinates to hg38
3. Extract regional LD matrix from GTEx genotypes
4. Intersect GWAS and eQTL SNPs
5. Run eCAVIAR colocalization test
6. Output results for each gene-tissue-GWAS-version combination

#### `countColoc.R`
**Results aggregation and filtering script**

**Functions:**
- Aggregates colocalization results across all genes for a tissue-cancer pair
- Applies significance filters:
  - CLPP (Colocalization Posterior Probability) ≥ 0.01
  - GWAS p-value ≤ 5×10⁻⁷
- Identifies top variants for each gene based on:
  - Maximum CLPP score
  - Minimum GWAS p-value

**Output:** Filtered summary tables with significant colocalization hits

#### `gather_eGenes.R`
**Data preparation script**

**Functions:**
- Scans tissue directories for eQTL results files
- Identifies expression genes (eGenes) with FDR < 0.05
- Merges with gene annotation data
- Creates input datasets for colocalization analysis

**Output:**
- `eGene_colocalization.RDS`: Gene metadata for colocalization
- `eGene_intersection_forplots.RDS`: eGene counts for visualization

## Workflow

```
1. runQTLTools.sh/R    → Perform eQTL mapping with QTLtools (48 jobs)
2. gather_eGenes.R     → Identify eGenes and prepare gene list
3. runColoc.sh/R       → Perform colocalization analysis (1536 jobs)
4. countColoc.sh/R     → Summarize and filter results (480 jobs)
```

## Data Sources

- **GTEx v8**: Expression and genotype data
- **Cancer GWAS**: 10 cancer types from various consortia
- **Reference**: hg19-to-hg38 liftover chain file
- **LD Reference**: GTEx WGS genotypes (838 samples)

## Software Dependencies

- **R packages**: data.table, rtracklayer, GenomicRanges, bigsnpr, tidyverse
- **External tools**: QTLtools, eCAVIAR, PLINK2
- **System**: LSF job scheduler

## Key Parameters

- **Colocalization window**: ±1Mb around gene boundaries
- **LD estimation**: GTEx WGS reference panel
- **Significance thresholds**:
  - eQTL FDR < 0.05
  - GWAS p < 5×10⁻8
  - CLPP ≥ 0.01

## Output Structure

```
/rsrch5/scratch/epi/bhattacharya_lab/GTEX_compare/
├── [tissue]/
│   └── coloc/
│       └── [cancer]/
│           └── [gene]_eQTL_coloc_[version]_[method].txt.gz
└── colocalization_results_[tissue]_[cancer]_[indices].tsv
```

## Usage Notes

- Pipeline designed for high-performance computing with LSF scheduler
- Handles multiple GTEx processing versions for robustness
- Implements coordinate liftover for GWAS compatibility
- Uses eCAVIAR method for Bayesian colocalization testing
- Results filtered for genome-wide significance and meaningful colocalization probability
