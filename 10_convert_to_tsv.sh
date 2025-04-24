#!/bin/bash

echo "[10 - VCF to TSV] Started at $(date)"

mkdir -p results_dna/tsv

for file in results_dna/filtered/filtered_*.vcf.gz; do
    base=$(basename $file .vcf.gz)
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t]\n' $file > results_dna/tsv/${base}.tsv
done

echo "[10 - VCF to TSV] Done at $(date)"
