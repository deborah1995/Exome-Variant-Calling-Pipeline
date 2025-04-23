#!/bin/bash
set -e

SAMPLE_LIST="samples.txt"
REF="../hg38_index/bwa/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUT_DIR="results_dna"
THREADS=4
LOG_DIR="logs"

mkdir -p $OUT_DIR $LOG_DIR

while read SAMPLE; do
  echo "[02 - BWA] Running for $SAMPLE at $(date)"
  mkdir -p $OUT_DIR/$SAMPLE

  bwa mem -t $THREADS $REF trimmed/${SAMPLE}_R1_trimmed.fastq.gz trimmed/${SAMPLE}_R2_trimmed.fastq.gz > $OUT_DIR/$SAMPLE/${SAMPLE}.sam

  echo "[02 - BWA] Done for $SAMPLE at $(date)"
done < $SAMPLE_LIST
