#!/bin/bash

# Define paths
REF_FASTA="/home/bastoni/hg38_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GFF3_FILE="/home/bastoni/working_directory/Homo_sapiens.GRCh38.113.gff3"
RESULTS_DIR="/home/bastoni/working_directory/results_dna"
TSV_DIR="${RESULTS_DIR}/tsv"

# Create TSV directory if not exists
mkdir -p "$TSV_DIR"

# Annotate SNPs
echo " Annotating SNPs with bcftools csq..."
bcftools csq -f "$REF_FASTA" \
  -g "$GFF3_FILE" \
  -p a \
  -Oz -o "${RESULTS_DIR}/annotated_with_csq_snps.vcf.gz" \
  "${RESULTS_DIR}/filtered_snps.vcf"

# Annotate INDELs
echo " Annotating INDELs with bcftools csq..."
bcftools csq -f "$REF_FASTA" \
  -g "$GFF3_FILE" \
  -p a \
  -Oz -o "${RESULTS_DIR}/annotated_with_csq_indels.vcf.gz" \
  "${RESULTS_DIR}/filtered_indels.vcf"

# Index both files
echo "🔧 Indexing VCFs..."
bcftools index "${RESULTS_DIR}/annotated_with_csq_snps.vcf.gz"
bcftools index "${RESULTS_DIR}/annotated_with_csq_indels.vcf.gz"

# Convert to TSV with CHROM, POS, REF, ALT, annotation (BCSQ), and genotypes
echo " Extracting SNPs to TSV..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT]\n' \
  "${RESULTS_DIR}/annotated_with_csq_snps.vcf.gz" > "${TSV_DIR}/annotated_snps_with_samples.tsv"

echo " Extracting INDELs to TSV..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT]\n' \
  "${RESULTS_DIR}/annotated_with_csq_indels.vcf.gz" > "${TSV_DIR}/annotated_indels_with_samples.tsv"

echo " All done! SNPs and INDELs annotated and converted to TSV format in $TSV_DIR"
