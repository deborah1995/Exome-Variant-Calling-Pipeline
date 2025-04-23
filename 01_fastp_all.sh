#!/bin/bash
set -e

# === Configuration ===
SAMPLE_LIST="samples.txt"
INPUT_DIR="../esoma"
OUT_DIR="trimmed"
LOG_DIR="logs"
THREADS=4

mkdir -p $OUT_DIR $LOG_DIR

while read SAMPLE; do
  echo "[01 - fastp] Running for $SAMPLE at $(date)"
  
  fastp \
    -i ${INPUT_DIR}/${SAMPLE}_R1_sub.fastq.gz \
    -I ${INPUT_DIR}/${SAMPLE}_R2_sub.fastq.gz \
    -o ${OUT_DIR}/${SAMPLE}_R1_trimmed.fastq.gz \
    -O ${OUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz \
    -h ${OUT_DIR}/${SAMPLE}_fastp.html \
    -j ${OUT_DIR}/${SAMPLE}_fastp.json \
    -w $THREADS

  echo "[01 - fastp] Done for $SAMPLE at $(date)"
done < $SAMPLE_LIST
