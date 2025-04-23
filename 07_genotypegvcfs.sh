#!/bin/bash
set -e

# === Config ===
DB_WORKSPACE="genomicsdb"
REF="../hg38_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUTPUT_VCF="joint_genotyped.vcf.gz"
LOGDIR="logs"
THREADS=4

mkdir -p $LOGDIR

echo "[07 - GenotypeGVCFs] Started at $(date)" | tee $LOGDIR/genotypegvcfs.log

gatk GenotypeGVCFs \
  -R $REF \
  -V gendb://$DB_WORKSPACE \
  -O $OUTPUT_VCF \
  --tmp-dir ./tmp \
  --native-pair-hmm-threads $THREADS \
  2>> $LOGDIR/genotypegvcfs.log

echo "[07 - GenotypeGVCFs] Done at $(date)" | tee -a $LOGDIR/genotypegvcfs.log
