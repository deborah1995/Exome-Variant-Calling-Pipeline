#!/bin/bash
set -e

# === Config ===
GVCF_DIR="results_dna"
SAMPLE_LIST="samples.txt"
BED="grch38_exons_merged.bed"
WORKSPACE="genomicsdb"
LOGDIR="logs"
THREADS=4

# Remove workspace if re-running
rm -rf $WORKSPACE
mkdir -p $WORKSPACE $LOGDIR

# === Prepare input list ===
INPUTS=""
while read SAMPLE; do
  INPUTS+=" -V ${GVCF_DIR}/${SAMPLE}/${SAMPLE}.g.vcf.gz"
done < $SAMPLE_LIST

# === Run GenomicsDBImport ===
echo "[06 - GenomicsDBImport] Started at $(date)" | tee $LOGDIR/genomicsdbimport.log

gatk --java-options "-Xmx4g" GenomicsDBImport \
  --genomicsdb-workspace-path genomicsdb_fixed \
  --intervals grch38_clean3.bed \
  -V results_dna/Sample1/Sample1.g.vcf.gz \
  -V results_dna/Sample2/Sample2.g.vcf.gz \
  -V results_dna/Sample3/Sample3.g.vcf.gz \
  -V results_dna/Sample4/Sample4.g.vcf.gz \
  --batch-size 50 \
  --reader-threads 4 \
  --merge-input-intervals \
  --disable-sequence-dictionary-validation

echo "[06 - GenomicsDBImport] Done at $(date)" | tee -a $LOGDIR/genomicsdbimport.log
