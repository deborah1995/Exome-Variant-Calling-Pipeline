#!/bin/bash
set -e

GTF="../hg38_index/Homo_sapiens.GRCh38.113.gtf"
BED_TMP="grch38_exons.bed"
BED_SORTED="grch38_exons_sorted.bed"
BED_FINAL="grch38_exons_merged.bed"

echo "[BED] Extracting exon regions from $GTF ..."
grep -P "\texon\t" $GTF | \
  awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5}' > $BED_TMP

echo "[BED] Sorting BED file ..."
sort -k1,1 -k2,2n $BED_TMP > $BED_SORTED

echo "[BED] Merging overlapping exons ..."
bedtools merge -i $BED_SORTED > $BED_FINAL

echo "[BED] Final BED saved as $BED_FINAL"
