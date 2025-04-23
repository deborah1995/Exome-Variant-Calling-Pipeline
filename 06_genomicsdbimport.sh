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

gatk GenomicsDBImport \
  --genomicsdb-workspace-path $WORKSPACE \
  --intervals $BED \
  $INPUTS \
  --batch-size 50 \
  --reader-threads $THREADS \
  2>> $LOGDIR/genomicsdbimport.log

echo "[06 - GenomicsDBImport] Done at $(date)" | tee -a $LOGDIR/genomicsdbimport.log
