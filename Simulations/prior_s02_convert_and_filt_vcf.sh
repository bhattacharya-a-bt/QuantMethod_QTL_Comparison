module load plink/2.00-alpha
module load tabix/0.2.6

PASS=1
DATA_DIR="/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/ldref/1KG"
OUT_DIR="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/1KG_vcf"
MAF_THRESH=0.05
SAMPLE_FILE="/rsrch5/scratch/epi/sthead/GTEx_gencode_comp/pass${PASS}/files_for_analysis/1kg_eur_500_sample_ids"

mkdir -p ${OUT_DIR}

cd ${OUT_DIR}

for chr in {1..22}
do
   plink2 --pfile ${DATA_DIR}/all_hg38 \
   --keep ${SAMPLE_FILE} \
   --chr ${chr} \
   --recode vcf \
   --maf ${MAF_THRESH} \
   --snps-only \
   --allow-extra-chr 0 \
   --out genos_1kg_eur_500_snps_maf_0.01_chr${chr}

done


# compress each vcf files

for chr in {1..22}
do
   bgzip genos_1kg_eur_500_snps_maf_0.01_chr${chr}.vcf

done


# index each compressed vcf file

for chr in {1..22}
do
   tabix -p vcf genos_1kg_eur_500_snps_maf_0.01_chr${chr}.vcf.gz

done

