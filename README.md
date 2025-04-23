# Exome-Variant-Calling-Pipeline
This repository contains all scripts and notes used to perform variant calling on a set of human exome sequencing samples (FASTQ). The goal is to identify genomic differences between samples derived from human liver biopsies. All analyses were conducted in a reproducible and modular manner.

Overview of the Analysis

The pipeline includes the following major steps:

Quality Control and Trimming using fastp

Alignment of reads to the human reference genome (GRCh38) using bwa mem

Sorting and Indexing using samtools

Marking Duplicates using Picard

Variant Calling per sample using GATK HaplotypeCaller

Joint Genotyping preparation using GATK GenomicsDBImport

🛠️ Tools and Environment

This pipeline was run on a Linux server with Conda. Tools used:

fastp

bwa

samtools

picard

gatk (v4.3.0.0)

Environment created with:

conda create -n exome_pipeline -c bioconda fastp bwa samtools gatk4 picard

⚙️ Script Descriptions

Each script processes all samples listed in samples.txt:

01_fastp_all.sh: quality filtering and adapter trimming.

02_bwa_align_all.sh: aligns reads to GRCh38 reference genome.

03_sort_index_all.sh: sorts and indexes BAM files.

04_mark_duplicates_all.sh: marks PCR duplicates.

05_haplotypecaller_all.sh: runs GATK HaplotypeCaller per sample in GVCF mode.

06_genomicsdbimport.sh: prepares for joint genotyping.

generate_exome_bed.sh: extracts merged exon intervals from a GTF file.

🔁 Reproducibility and Logging

Each step logs progress and output to a timestamped file in logs/.

Parallelization is handled with multi-threaded flags where supported.

Scripts follow a numerical naming convention to ensure ordered execution.

📌 Notes

The reference genome used is GRCh38 primary assembly.

The BED file for exon intervals was derived from GENCODE GTF.

Large intermediate files and output results are excluded via .gitignore.

📈 Next Steps

Joint genotyping with GenotypeGVCFs

Variant filtering and annotation (e.g. using VEP or ANNOVAR)

Exploratory analyses and visualization (PCA, heatmaps, etc.)
