#!/bin/bash
set -e

SAMPLE_LIST="samples.txt"
ROOTDIR="results_dna"
LOGDIR="logs"
REF="../hg38_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
THREADS=4

mkdir -p $LOGDIR

while read SAMPLE; do
  BAM_IN=$ROOTDIR/${SAMPLE}/${SAMPLE}.marked.bam
  GVCF_OUT=$ROOTDIR/${SAMPLE}/${SAMPLE}.g.vcf.gz
  LOG_FILE=$LOGDIR/${SAMPLE}_haplotypecaller.log

  echo "[05 - HaplotypeCaller] Calling variants for $SAMPLE at $(date)" | tee $LOG_FILE

  gatk HaplotypeCaller \
    -R $REF \
    -I $BAM_IN \
    -O $GVCF_OUT \
    -ERC GVCF \
    -bamout $ROOTDIR/${SAMPLE}/${SAMPLE}.bamout.bam \
    --native-pair-hmm-threads $THREADS 2>> $LOG_FILE

  echo "[05 - HaplotypeCaller] Done for $SAMPLE at $(date)" | tee -a $LOG_FILE

done < $SAMPLE_LIST
