# GTEx Data Processing Pipeline

This repository contains scripts for processing GTEx (Genotype-Tissue Expression) RNA-seq data through multiple annotation versions and analysis approaches. The pipeline supports processing data with GENCODE v27, v38, and v45 annotations using both Salmon and STAR alignment strategies.

## Overview

The pipeline processes GTEx BAM files through the following major steps:
1. **BAM to FASTQ conversion** - Extract paired-end reads from original BAM files
2. **Quantification/Alignment** - Process reads using either Salmon or STAR
3. **Gene-level summarization** - Aggregate transcript-level data to genes
4. **Normalization** - Apply TMM normalization for comparative analysis
5. **BED file generation** - Create QTLtools-compatible expression files

## Pipeline Components

### 1. BAM to FASTQ Conversion and Salmon Quantification

**Primary Script:** `bam_to_fastq_master.sh`
**Batch Submission:** `importGTExv45_all.sh`

This component:
- Converts GTEx BAM files to paired-end FASTQ format using samtools
- Runs Salmon quantification against GENCODE v45 transcriptome
- Processes samples in parallel using LSF job arrays (1-48 tissues)

**Key Features:**
- Name-sorts BAM files before FASTQ extraction
- Uses Salmon with bias correction and mapping validation
- Automatically cleans up intermediate files to save storage

### 2. STAR Alignment Pipeline

**Primary Scripts:** 
- `alignSTAR.sh` - Individual sample processing
- `alignSTAR_by_file.sh` - Batch processing by tissue file
- `alignSTAR_jobarray.sh` - LSF job array submission

This component:
- Performs STAR alignment against GRCh38 reference genome
- Generates coordinate-sorted BAM files with splice junction information
- Supports both GENCODE v38 and v45 processing

**STAR Parameters:**
- Multi-mapping: up to 20 alignments per read
- Splice junction overhang: minimum 8bp
- Intron size: 20bp to 1Mb
- Outputs coordinate-sorted BAM with standard attributes

### 3. Gene Counting with featureCounts

**Primary Script:** `countGene_featureCounts.R`

Features:
- Uses Rsubread::featureCounts for gene-level quantification
- Supports both GENCODE v38 and v45 annotations
- Processes STAR-aligned BAM files
- Handles paired-end reads with multi-threading

### 4. Transcript-Level Processing and Normalization

**Primary Script:** `importGTExv45.R`

This R script:
- Uses tximeta for transcript-level import from Salmon quantifications
- Summarizes transcript expression to gene-level using summarizeToGene()
- Applies TMM (Trimmed Mean of M-values) normalization using edgeR
- Creates multiple expression matrices:
  - Raw counts
  - TPM (Transcripts Per Million)
  - TMM-normalized CPM
  - Log2(TMM+1) transformed values

### 5. BED File Generation for QTL Analysis

**Primary Scripts:**
- `writeExpBed.R` - Main processing script
- `writeExpBed.sh` - LSF batch submission

Creates QTLtools-compatible BED files containing:
- Chromosome coordinates (chr1-22 only)
- Gene identifiers (Ensembl IDs)
- Log2-transformed expression values
- Proper BED format with bgzip compression and tabix indexing

**Processing Rules:**
- Filters genes with ≤0.1 TPM in >25% of samples
- Removes duplicate samples and genes
- Applies log2(x+1) transformation
- Generates separate files for v27, v38, and v45 annotations

## Key Dependencies

### R Packages
- **tximeta** - Transcript-level import and metadata
- **edgeR** - TMM normalization and differential expression
- **DESeq2** - Alternative normalization methods
- **SummarizedExperiment** - Data container objects
- **Rsubread** - Feature counting from alignments
- **rtracklayer** - GTF/GFF file handling
- **biomaRt** - Gene annotation and mapping

### External Tools
- **Salmon** - Transcript quantification
- **STAR** - Splice-aware alignment
- **samtools** - BAM file manipulation
- **bgzip/tabix** - File compression and indexing

## Usage Examples

### Process single tissue through Salmon pipeline:
```bash
# Submit job array for all tissues
bsub < importGTExv45_all.sh
```

### Generate BED files for QTL analysis:
```bash
# Process specific tissue (index 1-48)
bsub -J "writeExpBed[1-48]" writeExpBed.sh
```

### Run STAR alignment with gene counting:
```bash
# Submit combined STAR + featureCounts jobs
bsub < alignSTAR_jobarray.sh
```

## Output Files

### Expression Matrices
- **[TISSUE]_gene.RDS** - SummarizedExperiment with multiple assays
- **[TISSUE]_transcripts.RDS** - Transcript-level expression data

### QTL-ready Files
- **[TISSUE]_v[VERSION].bed.gz** - Compressed BED files with expression
- **[TISSUE]_v[VERSION].bed.gz.tbi** - Tabix indices

### Count Matrices
- **[TISSUE]_v[VERSION].RDS** - featureCounts output with annotations

## Notes

- The pipeline is designed for LSF batch systems with job arrays
- Memory requirements vary by tissue size (10-80GB typical)
- Storage requirements are substantial due to intermediate FASTQ files
- All scripts include automatic cleanup of temporary files
- Version compatibility is maintained across GENCODE releases v27, v38, and v45

## Quality Control

The pipeline includes several QC measures:
- Duplicate sample and gene removal
- Low-expression gene filtering (≤0.1 TPM in >25% samples)
- File size validation for alignment outputs
- Chromosome filtering (autosomes only: chr1-22)
- Proper handling of gene ID versioning across GENCODE releases
