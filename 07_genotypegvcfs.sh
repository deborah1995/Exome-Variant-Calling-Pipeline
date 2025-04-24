#!/bin/bash
set -euo pipefail

# === Config ===
REFERENCE="/home/bastoni/hg38_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENOMICSDB="genomicsdb_fixed"
OUTPUT_VCF="results_dna/genotyped.vcf.gz"
LOGFILE="logs/genotypegvcfs.log"

mkdir -p results_dna logs

# === Run GenotypeGVCFs ===
echo "[07 - GenotypeGVCFs] Started at $(date)" | tee $LOGFILE

gatk --java-options "-Xmx4g" GenotypeGVCFs \
  -R "$REFERENCE" \
  -V "gendb://$GENOMICSDB" \
  -O "$OUTPUT_VCF" \
  --disable-sequence-dictionary-validation \
  2>> $LOGFILE

echo "[07 - GenotypeGVCFs] Done at $(date)" | tee -a $LOGFILE

