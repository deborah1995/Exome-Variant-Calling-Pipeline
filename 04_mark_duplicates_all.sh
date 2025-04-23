#!/bin/bash
set -e

SAMPLE_LIST="samples.txt"
OUT_DIR="results_dna"
LOG_DIR="logs"

mkdir -p $LOG_DIR

while read SAMPLE; do
  echo "[04 - MarkDuplicates] Running for $SAMPLE at $(date)"

  picard AddOrReplaceReadGroups \
    I=$OUT_DIR/$SAMPLE/${SAMPLE}.sorted.bam \
    O=$OUT_DIR/$SAMPLE/${SAMPLE}.rg.bam \
    RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$SAMPLE

  picard MarkDuplicates \
    I=$OUT_DIR/$SAMPLE/${SAMPLE}.rg.bam \
    O=$OUT_DIR/$SAMPLE/${SAMPLE}.marked.bam \
    M=$OUT_DIR/$SAMPLE/${SAMPLE}.metrics.txt \
    CREATE_INDEX=true

  echo "[04 - MarkDuplicates] Done for $SAMPLE at $(date)"
done < $SAMPLE_LIST
