#!/bin/bash
set -e

SAMPLE_LIST="samples.txt"
OUT_DIR="results_dna"
LOG_DIR="logs"

mkdir -p $LOG_DIR

while read SAMPLE; do
  echo "[03 - Sort/Index] Running for $SAMPLE at $(date)"

  samtools view -bS $OUT_DIR/$SAMPLE/${SAMPLE}.sam | samtools sort -o $OUT_DIR/$SAMPLE/${SAMPLE}.sorted.bam
  samtools index $OUT_DIR/$SAMPLE/${SAMPLE}.sorted.bam

  echo "[03 - Sort/Index] Done for $SAMPLE at $(date)"
done < $SAMPLE_LIST
