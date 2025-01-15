# Mitochondrial-variants-detection-from-single-cell-RNA-sequencing-by-maegatk

* Purpose: Extract mitochondrial variants for cell lineage tracing
* Refer to https://github.com/petervangalen/MAESTER-2021, https://github.com/caleblareau/maegatk

inputs files:

* fastq
* scRNA-seq seurat
* cellranger count output

## Step 1: add cell barcodes and UMI from read 1 to read 2.

```bash
gunzip barcodes.tsv.gz #barcodes.tsv.gz from cellranger output

Rscript Assemble_fastq.R HPAP101/fastq HPAP101 barcodes.tsv 16 12 # Assemble_fastq.R can be found in the folder

fastqc HPAP101.fastq.gz # always run fastqc when a new fastq is generated

```

## Step 2: Trim 24bp from the start of Read 2 

```bash
homerTools trim -5 24 HPAP101.fastq.gz

fastqc HPAP101.fastq.gz.trimmed
```

## Step 3: Align to hg38 reference genome using STAR

```bash
# download the reference
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz

# generate idx

STAR --runMode genomeGenerate \
 --runThreadN 40  \
 --genomeDir /gpfs/research/scratch/lg23w/STAR/10X_human_reference/refdata-gex-GRCh38-2024-A/star2.7.9_idx \
 --genomeFastaFiles /gpfs/research/scratch/lg23w/STAR/10X_human_reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
 --sjdbGTFfile /gpfs/research/scratch/lg23w/STAR/10X_human_reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf

# star alignment
STAR --runMode alignReads  \
--runThreadN 40 \
--genomeDir /gpfs/research/scratch/lg23w/STAR/10X_human_reference/refdata-gex-GRCh38-2024-A/star2.7.9_idx \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn /gpfs/research/jwanggroup/jwang_group/FASTQ/sciRNAseq/HPAP101/pre-processing/HPAP101.fastq.gz.trimmed \
--outFileNamePrefix /gpfs/research/jwanggroup/jwang_group/FASTQ/sciRNAseq/HPAP101/pre-processing/trimmed_10X/HPAP101_trimmed_10X_

```

## Step 4: Add CB and UMI as sam tags

```bash
# generate bam index


```
