# GTEx TWAS Analysis Pipeline

This repository contains scripts for training transcriptome-wide association study (TWAS) models using GTEx data and performing TWAS analysis with cancer GWAS datasets.

## Overview

The pipeline builds gene expression prediction models from GTEx genotype and expression data, then uses these models to predict gene expression effects from GWAS summary statistics. This approach identifies genes whose predicted expression levels are associated with cancer risk, potentially revealing causal pathways from genetic variants to disease through gene expression.

## Files Description

### Shell Scripts

#### `trainTWAS.sh`
- **Purpose**: LSF job submission script for training TWAS prediction models
- **Job Array**: 1-576 jobs (48 tissues × 12 indices per tissue)
- **Resources**: 10 CPUs, 100GB memory, 24-hour time limit
- **Queue**: medium
- **Functions**:
  - Distributes model training across tissues and dataset versions
  - Each tissue gets 12 processing jobs for parallel execution
- **Calls**: `trainTWAS.R` with tissue name and secondary index

#### `runTWAS.sh`
- **Purpose**: LSF job submission script for performing TWAS analysis
- **Job Array**: 1-1000 jobs (100 model chunks × 10 cancer GWAS)
- **Resources**: 1 CPU, 40GB memory, 96-hour time limit
- **Queue**: long
- **Functions**:
  - Distributes TWAS testing across pre-trained models and GWAS datasets
  - Calculates array and GWAS indices from job index
- **Calls**: `runTWAS.R` with array and GWAS indices

### R Scripts

#### `trainTWAS.R`
**TWAS model training script**

**Key Functions:**
- Processes covariate-corrected gene expression data
- Trains prediction models using three methods:
  - **Elastic Net**: Regularized linear regression (α=0.5)
  - **BLUP**: Best Linear Unbiased Prediction
  - **SuSiE**: Sum of Single Effects fine-mapping
- Selects optimal model based on cross-validated R²
- Saves models with R² > 0.01 for downstream analysis

**Input Data:**
- **Expression**: Covariate-corrected BED files from QTLtools
- **Genotypes**: GTEx WGS reference panel (838 samples)
- **Versions**: v27, v38, v45 with STAR/Salmon processing

**Analysis Steps:**
1. Load and format expression data with strand information
2. Correct for covariates using QTLtools
3. Extract regional genotypes (±1Mb around genes)
4. Train multiple prediction models with 5-fold CV
5. Select best model based on prediction accuracy
6. Save model weights and metadata

**Output:** Individual RDS files containing model weights, R², and SNP annotations

#### `runTWAS.R`
**TWAS association testing script**

**Key Functions:**
- Loads pre-trained TWAS models and GWAS summary statistics
- Performs genomic coordinate liftover (hg19 → hg38)
- Computes TWAS Z-statistics using the formula:
  ```
  Z_TWAS = (W^T × Z_GWAS) / sqrt(W^T × LD × W)
  ```
  Where W = SNP weights, Z_GWAS = GWAS Z-scores, LD = linkage disequilibrium matrix

**Input Data:**
- **TWAS Models**: Pre-trained prediction models from `trainTWAS.R`
- **GWAS Files**: 10 cancer GWAS summary statistics (same as colocalization pipeline)
- **LD Reference**: GTEx genotype data for LD matrix computation

**Analysis Steps:**
1. Load GWAS data and perform coordinate liftover
2. Process TWAS models in batches (distributed across 100 array jobs)
3. Extract regional genotypes and compute LD matrix
4. Align alleles between GWAS and TWAS model
5. Calculate TWAS Z-statistic for gene-tissue-cancer combinations
6. Identify top GWAS signals in extended gene regions (±1Mb)

**Output:** TWAS association results with Z-scores and top regional GWAS p-values

## Workflow

```
1. trainTWAS.sh/R      → Train gene expression prediction models (576 jobs)
2. runTWAS.sh/R        → Perform TWAS association testing (1000 jobs)
```

## Data Sources

- **GTEx v8**: Expression and genotype data across 48 tissues
- **Cancer GWAS**: 10 cancer types (Endometrial, Renal, Prostate, Ovarian, Melanoma, Lung, Head & Neck, Glioma, Colorectal, Breast)
- **Processing Versions**: Multiple GTEx data processing pipelines (v27, v38, v45)
- **RNA-seq Methods**: STAR and Salmon alignment/quantification

## Software Dependencies

- **R packages**: 
  - `isotwas`: TWAS model training and testing
  - `data.table`, `bigsnpr`: Data manipulation and genotype handling
  - `rtracklayer`, `GenomicRanges`: Genomic coordinate operations
- **External tools**: QTLtools, PLINK2, bgzip, tabix
- **System**: LSF job scheduler

## Key Parameters

- **Training window**: ±1Mb around gene boundaries
- **Cross-validation**: 5-fold CV for model selection
- **R² threshold**: Minimum 0.01 for model retention
- **Elastic Net α**: 0.5 (equal L1/L2 penalty)
- **LD estimation**: GTEx WGS reference panel

## Model Information

Each TWAS model (RDS file) contains:
- **Feature**: Gene identifier (ENSG ID)
- **SNP**: Variant identifiers
- **Chromosome/Position**: Genomic coordinates (hg38)
- **Weight**: Prediction weights for each SNP
- **R²**: Cross-validated prediction accuracy
- **Alleles**: Reference and alternative alleles

## Usage Notes

- Pipeline requires substantial computational resources (hundreds of CPU-hours)
- Models are tissue and processing-method specific
- TWAS results provide gene-level association statistics
- Compatible with standard GWAS summary statistic formats
- Implements proper allele alignment and coordinate liftover
- Results can be integrated with colocalization analysis for comprehensive interpretation
