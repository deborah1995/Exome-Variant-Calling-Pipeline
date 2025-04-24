#!/bin/bash
set -euo pipefail

LOGDIR="logs"
mkdir -p $LOGDIR

echo "[08b - Select INDELs] Started at $(date)" | tee $LOGDIR/filter_indels.log

gatk SelectVariants \
  -R /home/bastoni/hg38_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -V results_dna/genotyped.vcf.gz \
  --select-type-to-include INDEL \
  -O results_dna/filtered_indels.vcf \
  2>> $LOGDIR/filter_indels.log

echo "[08b - VariantFiltration INDELs] Filtering started at $(date)" | tee -a $LOGDIR/filter_indels.log

gatk VariantFiltration \
  -R /home/bastoni/hg38_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -V results_dna/filtered_indels.vcf \
  -O results_dna/filtered/filtered_indels.vcf \
  --filter-expression "QD < 2.0 || FS > 200.0" \
  --filter-name "basic_indel_filter" \
  2>> $LOGDIR/filter_indels.log

bgzip -f results_dna/filtered/filtered_indels.vcf
tabix -p vcf results_dna/filtered/filtered_indels.vcf.gz

echo "[08b - INDELs filtered] Done at $(date)" | tee -a $LOGDIR/filter_indels.log
