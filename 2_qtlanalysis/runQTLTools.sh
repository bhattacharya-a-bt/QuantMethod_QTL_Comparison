#!/bin/bash
#BSUB -J "qtltools[1-48]"
#BSUB -o /rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/logs/qtltools_%J_%I.out
#BSUB -e /rsrch5/home/epi/abhattacharya3/processGTEx/2_qtlanalysis/logs/qtltools_%J_%I.err
#BSUB -q long
#BSUB -W 200:00
#BSUB -n 1
#BSUB -M 80
#BSUB -R "span[hosts=1]"

# Create log directory if it doesn't exist
mkdir -p /rsrch5/scratch/epi/abhattacharya3/GTEx_compare/logs

# Load QTLtools module
module load qtltools

# Get the ith folder from the directory
BASE_DIR="/rsrch5/scratch/epi/bhattacharya_lab"
FOLDERS=($(ls -d ${BASE_DIR}/*/))
TISSUE_DIR="${FOLDERS[$((LSB_JOBINDEX - 1))]}"
TISSUE=$(basename "$TISSUE_DIR")

echo "Processing job ${LSB_JOBINDEX}: TISSUE_DIR=${TISSUE_DIR}, TISSUE=${TISSUE}"

# Set file paths
VCFFILE="/rsrch5/scratch/epi/abhattacharya3/GTEx_compare/GTEx.WGS.838.vcf.gz"
COVARFILE=$(find "$TISSUE_DIR" -name "*_formatted_covariates.txt" | head -1)

# Check if covariate file was found
if [[ -z "$COVARFILE" ]]; then
    echo "ERROR: No covariate file matching *_formatted_covariates.txt found in $TISSUE_DIR"
    exit 1
fi

echo "VCF file: $VCFFILE"
echo "Covariate file: $COVARFILE"

# Check if input files exist
if [[ ! -f "$VCFFILE" ]]; then
    echo "ERROR: VCF file not found: $VCFFILE"
    exit 1
fi

if [[ ! -f "$COVARFILE" ]]; then
    echo "ERROR: Covariate file not found: $COVARFILE"
    exit 1
fi

# Find all bed.gz files in the tissue directory
BED_FILES=($(find "$TISSUE_DIR" -name "*.bed.gz"))

if [[ ${#BED_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No .bed.gz files found in $TISSUE_DIR"
    exit 1
fi

echo "Found ${#BED_FILES[@]} BED files in $TISSUE_DIR"

# Process each BED file
for BEDFILE in "${BED_FILES[@]}"; do
    echo "Processing BED file: $BEDFILE"
    
    # Extract prefix by removing .bed.gz extension
    BASENAME=$(basename "$BEDFILE")
    PREFIX="${BASENAME%.bed.gz}"
    
    # Set output files
    OUTFILE="${TISSUE_DIR}/${PREFIX}_cisQTL.txt"
    OUTFILE_NOM="${TISSUE_DIR}/${PREFIX}_cisQTL_nominal.txt"
    
    echo "  PREFIX: $PREFIX"
    echo "  Output file: $OUTFILE"
    echo "  Nominal output file: $OUTFILE_NOM"
    
    # Generate random seed (using job index, current time, and file hash for uniqueness)
    SEED=$((LSB_JOBINDEX * $(date +%s) * $(stat -c%s "$BEDFILE") % 2147483647))
    echo "  Using SEED: $SEED"
    
    # Run QTLtools with permutations
    echo "  Starting QTLtools permutation analysis for $PREFIX..."
    QTLtools cis --normal \
        --vcf "$VCFFILE" \
        --bed "$BEDFILE" \
        --cov "$COVARFILE" \
        --permute 1000 \
        --seed "$SEED" \
        --out "$OUTFILE"
    
    if [[ $? -ne 0 ]]; then
        echo "ERROR: QTLtools permutation analysis failed for $PREFIX"
        exit 1
    fi
    
    # Run QTLtools nominal analysis
    echo "  Starting QTLtools nominal analysis for $PREFIX..."
    QTLtools cis --normal \
        --vcf "$VCFFILE" \
        --bed "$BEDFILE" \
        --cov "$COVARFILE" \
        --nominal 1 \
	--std-err \
        --out "$OUTFILE_NOM"
    
    if [[ $? -ne 0 ]]; then
        echo "ERROR: QTLtools nominal analysis failed for $PREFIX"
        exit 1
    fi
    
    echo "  QTLtools analysis completed successfully for $PREFIX"
done

echo "Job ${LSB_JOBINDEX} completed successfully for tissue $TISSUE"