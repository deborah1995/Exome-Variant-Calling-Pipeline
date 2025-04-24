#!/bin/bash
set -euo pipefail

LOGDIR="logs"
mkdir -p $LOGDIR

echo "[08a - Select SNPs] Started at $(date)" | tee $LOGDIR/filter_snps.log

gatk SelectVariants \
  -R /home/bastoni/hg38_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -V results_dna/genotyped.vcf.gz \
  --select-type-to-include SNP \
  -O results_dna/filtered_snps.vcf \
  2>> $LOGDIR/filter_snps.log

echo "[08a - VariantFiltration SNPs] Filtering started at $(date)" | tee -a $LOGDIR/filter_snps.log

gatk VariantFiltration \
  -R /home/bastoni/hg38_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -V results_dna/filtered_snps.vcf \
  -O results_dna/filtered/filtered_snps.vcf \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "basic_snp_filter" \
  2>> $LOGDIR/filter_snps.log

bgzip -f results_dna/filtered/filtered_snps.vcf
tabix -p vcf results_dna/filtered/filtered_snps.vcf.gz

echo "[08a - SNPs filtered] Done at $(date)" | tee -a $LOGDIR/filter_snps.log
