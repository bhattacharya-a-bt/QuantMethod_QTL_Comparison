module load plink/2.00-alpha
module load tabix/0.2.6

GENO_PASS=1
DATA_FILE="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG/all_hg38" # path to plink pfile
OUT_DIR="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${GENO_PASS}/files_for_analysis/1KG_vcf" # output directory
MAF_THRESH=0.05 # minor allele frequency threshold
SAMPLE_FILE="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${GENO_PASS}/files_for_analysis/1kg_eur_500_sample_ids" # output from prior_s01
OUT_FILE="genos_1kg_eur_500_snps_maf_0.01" # basename of output file

mkdir -p ${OUT_DIR}

cd ${OUT_DIR}

for chr in {1..22}
do
   plink2 --pfile ${DATA_FILE} \
   --keep ${SAMPLE_FILE} \
   --chr ${chr} \
   --recode vcf \
   --maf ${MAF_THRESH} \
   --snps-only \
   --allow-extra-chr 0 \
   --out ${OUT_FILE}_chr${chr}

done


# compress each vcf files

for chr in {1..22}
do
   bgzip ${OUT_FILE}_chr${chr}.vcf

done


# index each compressed vcf file

for chr in {1..22}
do
   tabix -p vcf ${OUT_FILE}_chr${chr}.vcf.gz

done

