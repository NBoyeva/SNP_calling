# Bash Pipeline for SNP Calling and Annotation

## General Info
A pipeline for SNP analysis. Includes FASTQ automatic preprocessing, SNP calling, and annotation.  
First, FASTQ files are trimmed based on quality control with fastp. Then alignments on a reference genome are made separately for each lane of the sample. Reference genome is indexed if necessary. Then alignments are sorted, duplicates are removed, read groups are assigned.  
Alignments are merged for each sample, merged file is indexed, and SNP calling is performed.  
VCF files retrieved are annotated, filtered by quality, and optionally filtered by ClinVar database.   
Additional TSV report is retrieved for specific columns.

## Prerequisites
* samtools (v.1.13)
* fastp (v.0.23.2)
* Bowtie2 (v.2.4.4)
* GATK (v.4.5.0.0)
* SnpEff & SnpSift
* Python (v.3.11.5)
* Java (v.17.0.10)

## Launch
### Arguments

| Flag | Description |
| --- | --- |
| -i | Path to directiry with input FASTQ files |
| -o | Path to output directory  |
| -n | Number of cores to be involved during pipeline operating time |
| -d | Path to the directory with FASTA reference file  |
| -f | Reference FASTA filename  |
| -b | Basename of the indexes for reference FASTA  |
| -g | Path to the GATK .py file  |
| -s | Path to the folder with SnpEff and SnpSift .jar file  |

### Example of Usage
./pipeline.sh -i data \  
-o output \  
-d reference \  
-f human_genome.fasta \  
-b human_genome \  
-g /home/user/gatk-4.5.0.0/gatk \  
-s /home/user/soft/snpEff \  
-r 11

## Workflow Scheme
![snp scheme](https://github.com/NBoyeva/SNP_calling/assets/149397882/b6fda9e3-9fcf-44b6-aa8a-4f7e5b7ab3cd)
