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
samtools index HPAP101_trimmed_10X_Aligned.sortedByCoord.out.bam

# Add CB and UMI as sam tag with Tag_CB_UMI.sh

INPUT=HPAP101_trimmed_10X_Aligned.sortedByCoord.out.bam
OUTPUT="$(echo "${INPUT/.bam/.10x.bam}")"
/gpfs/home/lg23w/mitochondrial/Tag_CB_UMI.sh $INPUT  $OUTPUT
echo "Job completed!"

# output is HPAP101_trimmed_10X_Aligned.sortedByCoord.out.10x.bam
```
## Step 5: Subset for chrM, and filtering CBs, which based on
* Good quality CBs from scRNAseq analysis by ncount, nfeature
* contain >100 MT reads 

```bash
# subset for chrM
samtools view HPAP101_trimmed_10X_Aligned.sortedByCoord.out.10x.bam chrM -b > HPAP101_10X_MT.bam

# filtering, 1. export the CBs from Seurat as HPAP101_scRNA_barcode.txt (after filtering by ncount, nfeature) 2. CBs >100 MT reads, to get that,

## Count the reads map to each CB
samtools view HPAP101_10X_MT.bam | grep -o 'CB:Z:[^[:space:]]*' | sort | uniq -c | sort -nr > HPAP101_10X_MT_CB_output.txt

## Create a filtered list of CB
awk '$1 >= 100 {print $2}' HPAP101_10X_MT_CB_output.txt > MT_filtered_barcodes.txt

### test if the format is correct before proceeding
cat -A MT_filtered_barcodes.txt | head -n 20
wc -l MT_filtered_barcodes.txt ##### for HPAP101, the output is 2743
samtools view HPAP101_10X_MT.bam | grep "$(head -n 1 MT_filtered_barcodes.txt)" 

## to get the common set between MT_filtered_barcodes.txt and HPAP101_scRNA_barcode.txt. To do that, we need to add CB:Z: to HPAP101_scRNA_barcode.txt to make same format.

wc -l HPAP101_scRNA_barcode.txt #####3252 for HPAP101
awk '{print "CB:Z:"$1}' HPAP101_scRNA_barcode.txt > HPAP101_scRNA_barcode_match.txt

## Now we need to get the common set between MT CBs and scRNAseq CBs
grep -Fwf pre-processing/trimmed_10X/HPAP101_10X_MT_CB_output.txt cellranger_output/HPAP101_scRNA_barcode_match.txt > common_CB.txt
wc -l common_CB.txt  ###### for HPAP101, the output is 2716.

## Finally, we can use the common CBs to filter the bam file
samtools view -H HPAP101_10X_MT.bam > headers.sam
samtools view HPAP101_10X_MT.bam | grep -Ff common_CB.txt | cat headers.sam - | samtools view -b -o filtered_HPAP101_10X_MT_common.bam

### sorting, mark duplicated reads and then removed
samtools sort -o sorted_filtered_HPAP101_10X_MT_common.bam filtered_HPAP101_10X_MT_common.bam
samtools markdup -r sorted_filtered_HPAP101_10X_MT_common.bam sorted_filtered_HPAP101_10X_MT_common_dedup.bam

### change the header for further merge
samtools view -H sorted_filtered_HPAP101_10X_MT_common_dedup.bam | sed -e 's/SN:chr/SN:/g' -e 's/SN:M/SN:MT/' > header_modified.sam
samtools reheader header_modified.sam sorted_filtered_HPAP101_10X_MT_common_dedup.bam > renamed_sorted_filtered_HPAP101_10X_MT_common_dedup.bam


### for scRNA_seq bam
samtools view -H possorted_genome_bam.bam > headers.sam ## possorted_genome_bam.bam is from 10X genomics
samtools view possorted_genome_bam.bam | grep -Ff common_CB.txt | cat headers.sam - | samtools view -b -o filtered_HPAP101_10X_scRNA_common.bam

## merge (not nessary if only exam the mitochondrial results)
## Note, the headers of these two bam files are different, make sure the headers are matched before merge
samtools merge HPAP101_common.bam filtered_HPAP101_10X_scRNA_common.bam ../pre-processing/trimmed_10X/filtered_HPAP101_10X_MT_common_1.bam

```

## Step 6: Run `maegatk` to ge the .rds output
```bash
maegatk bcall -i /gpfs/research/jwanggroup/jwang_group/FASTQ/sciRNAseq/HPAP101/pre-processing/trimmed_10X/sorted_filtered_HPAP101_10X_MT_common_dedup.bam -o maester_dedup_qc -z -qc
```

# Downstream analysis 

## coverage
```r
# Load data
HPAP101_maegatk=readRDS("/Users/liguo/Desktop/mitochondria/HPAP101/maegatk_HPAP101.rds")
HPAP101_seurat=readRDS("/Users/liguo/Desktop/mitochondria/HPAP101/HPAP101_040424.rds")
HPAP101_seurat_common=HPAP101_seurat[,colnames(HPAP101_maegatk)]

# Set y axis parameters
ymax <- 100

# Gene locations
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
  mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(ymax*1.2,ymax*1.1), length.out = 15))

#Plot coverage Per base
base.tib <- tibble(base = 1:16569, depth = rowMeans(assays(HPAP101_maegatk)[["coverage"]]))

ggplot(base.tib) +
  geom_bar(aes(x = base, y = ifelse(depth > 1, yes = depth, no = NA)),
           stat = "identity", fill = "#64b53b", width = 1) +
  coord_cartesian(ylim = c(1, ymax), xlim = c(700, 15900)) +
  scale_y_continuous(trans = "log10") +
  geom_segment(data = GenePos.tib, aes(x = start, y = ycoord, xend = end, yend = ycoord)) +
  geom_text(data = GenePos.tib, aes(x = mid, y = ycoord - ymax * 0.2, label = Names), size = 3) +
  ylab("Mean coverage per cell") + xlab("Position along chrM") +
  theme_classic() +
  theme(aspect.ratio = 0.5)

# Plot Mean Coverage For Top Cells
cells.tib <- tibble(
  cell = colnames(HPAP101_maegatk),
  depth = HPAP101_maegatk$depth
)
topcells.tib <- cells.tib %>% slice_max(order_by = depth, n = 500)

ggplot(topcells.tib, aes(x = 1, y = depth)) +
  geom_violin() +
  geom_sina(size = 0.3) +
  coord_cartesian(ylim = c(0.1, 800)) +
  scale_y_continuous(trans = "log10") +
  ylab("Mean coverage per cell") + xlab("") +
  annotate("text", x = 1, y = max(topcells.tib$depth) * 1.5,
           label = round(mean(topcells.tib$depth), 2)) +
  theme_classic() +
  theme(aspect.ratio = 2, plot.title = element_text(hjust = 0.5)) +
  ggtitle("Mean coverage of top 500 cells")

# Find and Plot the top 5000 bases based on coverage depth
top.tib <- base.tib %>% arrange(desc(depth)) %>% mutate(key = row_number(), .before = 1) %>% filter(key %in% 1:5000)

ggplot(top.tib) +
  geom_bar(aes(x = key, y = depth), stat = "identity", fill = "#64b53b", width = 1) +
  coord_cartesian(ylim = c(1, ymax)) +
  scale_y_continuous(trans = "log10") +
  geom_label(aes(x = 2500, y = mean(top.tib$depth), label = round(mean(top.tib$depth), 2)),
             fill = "#64b53b") +
  ylab("Mean coverage per cell") + xlab("Rank sorted position") +
  ggtitle("Mean coverage of top 5000 bases") +
  theme_classic() +
  theme(aspect.ratio = 2)

DimPlot(HPAP101_seurat_common,reduction ="umap")

DimPlot(HPAP101_seurat_common,reduction ="umap",group.by = c("celltype","barcode"))
```

## cell_variant_calling
```r
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(tidyverse)
library(SummarizedExperiment)
library(Seurat)
library(data.table)
library(Matrix)
library(ComplexHeatmap)
library(readxl)
library(ggrastr)

install.packages("Cairo", type = "source")

# Load the data
HPAP101_maegatk=readRDS("/Users/liguo/Desktop/mitochondria/HPAP101/maegatk_HPAP101.rds")
HPAP101_seurat=readRDS("/Users/liguo/Desktop/mitochondria/HPAP101/HPAP101_040424.rds")
HPAP101_seurat_common=HPAP101_seurat[,colnames(HPAP101_maegatk)]

# Prepare allele frequency matrix. 
# Rows represents a position along the mitochondrial genome and the three possible disagreements with the reference (except 3107 has four possible disagreements because the reference is N) 

source("/Users/liguo/Desktop/mitochondria/HPAP101/MAESTER-2021/Auxiliary_files/210215_FunctionsGeneral.R")
af.dm <- data.matrix(computeAFMutMatrix(HPAP101_maegatk))*100
HPAP101_seurat_common$orig.ident=rownames(HPAP101_seurat_common@meta.data)

# Collect cell information
cells.tib <- tibble(cell = colnames(HPAP101_maegatk),
                    orig.ident = HPAP101_seurat_common$barcode,
                    CellType_RNA = HPAP101_seurat_common$celltype,
                    UMAP_1 = HPAP101_seurat_common@reductions$umap@cell.embeddings[,"UMAP_1"],
                    UMAP_2 = HPAP101_seurat_common@reductions$umap@cell.embeddings[,"UMAP_2"],
                    Mean_Cov = HPAP101_maegatk$depth)

# Make groups of Cell cell types
Cell_celltype.ls <- list(unionCells = cells.tib$cell,
                       alpha = filter(cells.tib, CellType_RNA == "alpha")$cell,
                       beta = filter(cells.tib, CellType_RNA == "beta")$cell,
                       delta = filter(cells.tib, CellType_RNA == "delta")$cell,
                       pp = filter(cells.tib, CellType_RNA == "pp")$cell,
                       exocrine = filter(cells.tib, CellType_RNA == "exocrine")$cell,
                       fibroblast = filter(cells.tib, CellType_RNA == "fibroblast")$cell,
                       endothelial = filter(cells.tib, CellType_RNA == "endothelial")$cell,
                       immune = filter(cells.tib, CellType_RNA == "immune")$cell)

# Get the mean allele frequency and coverage for every cell subset
mean_af.ls <- lapply(Cell_celltype.ls, function(x) rowMeans(af.dm[,x]))
mean_cov.ls <- lapply(Cell_celltype.ls, function(x) rowMeans(assays(HPAP101_maegatk)[["coverage"]][,x])[as.numeric(cutf(rownames(af.dm), d = "_"))])
names(mean_af.ls) <- paste0("mean_af.", names(mean_af.ls))
names(mean_cov.ls) <- paste0("mean_cov.", names(mean_cov.ls))

# Get the quantiles of the VAFs of each variant in each cell subset
quantiles <- c("q01" = 0.01, "q10" = 0.1, "q50" = 0.5, "q90" = 0.9, "q99" = 0.99)
start_time <- Sys.time()
quantiles.ls <- lapply(quantiles, function(x) lapply(Cell_celltype.ls, function(y) apply(af.dm[,y], 1, quantile, x) ))
Sys.time() - start_time

# Get the mean quality for each variant.
assays.ls <- lapply(HPAP101_maegatk@assays@data, function(x) as.matrix(x))
start_time <- Sys.time()
qual.num <- sapply(rownames(af.dm), function(x) {
  #x <- "2141_T>C"
  pos <- as.numeric( cutf(x, d = "_") )
  mut <- cutf(x, d = ">", f = 2)
  # Get the mean quality of reads for this call (only use cells in which the base was sequenced) - forward
  covered_fw <- assays.ls[[str_c(mut, "_counts_fw")]][pos,] > 0
  qual_fw <- assays.ls[[str_c(mut, "_qual_fw")]][pos, covered_fw]
  # Same for reverse
  covered_rev <- assays.ls[[str_c(mut, "_counts_rev")]][pos,] > 0
  qual_rev <- assays.ls[[str_c(mut, "_qual_rev")]][pos, covered_rev]
  qual <- mean(c(qual_fw, qual_rev))
  return(qual)
})
Sys.time() - start_time

# Collect all information in a tibble
vars.tib <- as_tibble(do.call(cbind, c(mean_af.ls, mean_cov.ls, unlist(quantiles.ls, recursive = F))), rownames = "var")
vars.tib <- add_column(vars.tib, quality = qual.num, .before = 2)

# Save for fast loading next time
write_tsv(vars.tib, "HPAP101_celltype_Variants1.txt")
vars.tib <- read_tsv("HPAP101_celltype_Variants1.txt")

# Variants 
voi.ch <- vars.tib %>% filter(mean_cov.unionCells > 20,
                              quality > 30) %>% .$var
voi.ch_1 <- vars.tib %>% filter(quality > 30) %>% .$var

# Assess transitions vs. transversions
str_view(voi.ch, "G>A|A>G|C>T|T>C"); mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )
vars.tib %>%
  ggplot(aes(x = mean_af.alpha, y = mean_af.beta, color = var %in% voi.ch)) +
  geom_point(shape = 1, size = 2)+theme_classic()

# Generate matrices with coverage, allele frequency and reference / variant reads
cov_voi.mat <- assays(HPAP101_maegatk)[["coverage"]][as.numeric(cutf(voi.ch, d = "_")),]
af_voi.mat <- af.dm[voi.ch,]

### Plot VAF in each cell, sorted from low to high, illustrating the selection process
popcol.df <- read_excel("/Users/liguo/Desktop/mitochondria/HPAP101/MAESTER-2021/Auxiliary_files/MAESTER_colors.xlsx")
pdf(paste0("HPAP101", "_1_Sorted_VAFs.pdf"), width = 10, height = 5)
par(mar=c(4,4,2,8), xpd = T)
plot(NA, xlim = c(0, ncol(af_voi.mat)), ylim = c(0, 100), xlab = "Cells (sorted separately for each variant)", ylab = "Variant Allele Frequency")
for (n in 1:length(voi.ch)) {
  v <- voi.ch[n]
  points(as.numeric(sort(af_voi.mat[v,])), pch = 16, col = popcol.df$hex[n])
  text(x = 1.1*ncol(af_voi.mat), y = 110-10*n, label = v, col = popcol.df$hex[n], pos = 4)
}
dev.off()

# Add coverage and allele frequency info from variants of interest to cells.tib
for (voi in voi.ch) {
  cells.tib <- cells.tib %>%
    left_join(as_tibble(assays(HPAP101_maegatk)[["coverage"]][as.numeric(cutf(voi, d = "_")),], rownames = "cell"), by = "cell") %>%
    left_join(as_tibble(af.dm[voi,], rownames = "cell"), by = "cell") %>%
    rename(value.x = str_c("cov_", str_replace(voi, ">", ".")), value.y = str_c("af_", str_replace(voi, ">", ".")))
}

### For each variant of interest, plot UMAP of cells, colored by VAF
heatcol.ch <- read_excel("/Users/liguo/Desktop/mitochondria/HPAP101/MAESTER-2021/Auxiliary_files/MAESTER_colors.xlsx", sheet = 2, col_names = "heatcol")$heatcol

pdf(paste0("HPAP101", "_2_VAF_UMAP.pdf"))
for (voi in voi.ch) {
  #voi <- voi.ch[1]
  cov_colname <- str_c("cov_", str_replace(voi, ">", "."))
  af_colname <- str_c("af_", str_replace(voi, ">", "."))
  
  # Select cells with three genotyped transcripts, then plot
  print(
    cells.tib %>% filter(.[[cov_colname]] > 3) %>%
      ggplot(aes_string(x = "UMAP_1", y = "UMAP_2", color = af_colname)) +
      geom_point_rast() +
      scale_color_gradientn(colors = heatcol.ch[2:10]) +
      theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
      ggtitle(voi)
  )
}
dev.off()

for (voi in voi.ch) {
  #voi <- voi.ch[1]
  cov_colname <- str_c("cov_", str_replace(voi, ">", "."))
  af_colname <- str_c("af_", str_replace(voi, ">", "."))
  
  # Select cells with three genotyped transcripts, then plot
  print(
    cells.tib %>% filter(.[[cov_colname]] > 3) %>%
      ggplot(aes_string(x = "UMAP_1", y = "UMAP_2", color = af_colname)) +
      geom_point_rast() +
      scale_color_gradientn(colors = heatcol.ch[2:10]) +
      theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
      ggtitle(voi)
  )
}

plot_data_list <- list()

for (voi in voi.ch) {
  cov_colname <- str_c("cov_", str_replace(voi, ">", "."))
  af_colname <- str_c("af_", str_replace(voi, ">", "."))
  
  # Filter cells and add the variant name
  plot_data_list[[voi]] <- cells.tib %>%
    filter(.[[cov_colname]] > 3) %>%
    mutate(VAF = .[[af_colname]], Variant = voi) %>%  # Add VAF and Variant columns
    select(UMAP_1, UMAP_2, VAF, Variant)
}

# Combine all variant data into one tibble
plot_data <- bind_rows(plot_data_list)

ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = VAF)) +
  geom_point(size = 0.5) +  # Adjust point size for clarity
  scale_color_gradientn(colors = heatcol.ch[2:10]) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 6),  # Adjust facet label size
    axis.text = element_blank(),         # Optional: Remove axis text
    axis.ticks = element_blank()         # Optional: Remove axis ticks
  ) +
  facet_wrap(~ Variant, nrow = 6, ncol = 20) +
  ggtitle("UMAP Plots for 117 Variants")


####### Generate 12 VAF per plot, totally 10
variant_groups <- split(voi.ch, ceiling(seq_along(voi.ch) / 10))

# Loop through each group of variants
for (i in seq_along(variant_groups)) {
  group <- variant_groups[[i]]
  
  # Combine data for the current group
  group_data <- bind_rows(
    lapply(group, function(voi) {
      cov_colname <- str_c("cov_", str_replace(voi, ">", "."))
      af_colname <- str_c("af_", str_replace(voi, ">", "."))
      
      cells.tib %>%
        filter(.[[cov_colname]] > 3) %>%
        mutate(VAF = .[[af_colname]], Variant = voi) %>%
        select(UMAP_1, UMAP_2, VAF, Variant)
    })
  )
  
  # Generate the plot for the current group
  p <- ggplot(group_data, aes(x = UMAP_1, y = UMAP_2, color = VAF)) +
    geom_point(size = 0.5,alpha=0.5) +
    scale_color_gradientn(colors = heatcol.ch) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      strip.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    ) +
    facet_wrap(~ Variant, nrow = 2, ncol = 5) +
    ggtitle(paste("UMAP Plots for Variants (Page", i, ")"))
  
  # Save the plot as a PDF
  ggsave(paste0("UMAP_VAF_Group_", i, ".pdf"), plot = p, width = 10, height = 4)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### heatmap ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Select qualify cells
cells.ch <- cells.tib %>% 
  filter(orig.ident != "undetermined") %>%
  select(cell,orig.ident,starts_with("cov_")) %>%
  filter(apply(.[,-c(1,2)], 1, function(x) all(x > 3))) %>%
  group_by(orig.ident) %>%
  filter(n()>3) %>%
  ungroup() %>%
  .$cell
 
table(cells.tib %>% filter(cell %in% cells.ch) %>% pull(orig.ident))

Islets=unique(cells.tib %>% filter(cell %in% cells.ch) %>% pull(orig.ident))
set.seed(123)
islet_colors <- setNames(colors()[sample(length(colors()), length(Islets))], Islets)
celliselts.ha=HeatmapAnnotation(CellIslet=cells.tib$orig.ident[cells.ch],
                                col=list(CellIslet=islet_colors))


hm1 <- pheatmap(af_voi.mat[,cells.ch],
                clustering_method = "ward.D2",
                color=c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D"),
                show_colnames = F, show_rownames = T)

af_values <- as.vector(af_voi.mat[, cells.ch])
breaks <- quantile(af_values, probs = seq(0, 1, length.out = 101), na.rm = TRUE)
af_percentile_mapped <- apply(af_voi.mat[, cells.ch], c(1, 2), function(x) {
  findInterval(x, vec = breaks, rightmost.closed = TRUE)
})
color_gradient <- colorRampPalette(c("#FCFCFC", "#FFDF5F", "#F14C2B", "#A31D1D"))(100)

pheatmap(
  af_percentile_mapped,
  clustering_method = "ward.D2",
  color = color_gradient,
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Percentile-Based Heatmap"
)

hm2 <- Heatmap(
  af_percentile_mapped,
  clustering_method = "ward.D2",
  color = color_gradient,
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Percentile-Based Heatmap",
  top_annotation = celltype.ha,
  border=T
)

hm2 <- pheatmap(
  af_percentile_mapped,
  color = color_gradient,
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Percentile-Based Heatmap",
  top_annotation = celltype.ha
)

hm3 <- pheatmap(
  af_percentile_mapped,
  color = color_gradient,
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Percentile-Based Heatmap",
  top_annotation = celliselts.ha
)

celltype_levels <- c("alpha", "beta", "delta", "pp", "exocrine", "fibroblast", "endothelial", "immune")
cells.tib_1=subset(cells.tib,cell %in% cells.ch)
cells.tib_1$CellType_RNA <- factor(cells.tib_1$CellType_RNA, levels = celltype_levels)

cells.tib_1 <- cells.tib_1[order(cells.tib_1$CellType_RNA), ]
cell_order <- cells.tib_1$cell

celltype.ha=HeatmapAnnotation(CellType_RNA=cells.tib_1$CellType_RNA,
                              col=list(CellType_RNA=c("alpha"="#4682B4",
                                                      "beta"="#CD5C5C", 
                                                      "delta"="#5F9EA0", 
                                                      "pp"="#87CEEB",
                                                      "exocrine"= "#FF8C00",
                                                      "fibroblast"="#FFD700",
                                                      "endothelial"="#7B68EE",
                                                      "immune"="#FF6347")))



hm4 <- pheatmap(
  af_voi.mat[,cells.tib_1$cell],
  color = color_gradient,
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "AF Heatmap",
  cluster_cols = FALSE,
  top_annotation = celltype.ha
)

Heatmap(
  af_percentile_mapped[, cells.ch[celltype_order]], # Reordered columns
  col = color_gradient,
  clustering_method = "ward.D2", 
  show_column_names = FALSE,
  show_row_names = TRUE,
  top_annotation = celltype.ha,
  border = TRUE,
  name = "Percentile-Based Heatmap"
)

hm5 <- pheatmap(
  af_percentile_mapped[,cells.tib_1$cell],
  color = color_gradient,
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Percentile-Based Heatmap",
  top_annotation = celltype.ha,
  cluster_cols = FALSE
)

cells.tib_2 <- cells.tib_1 %>%
  arrange(orig.ident)

head(cells.tib_2)

celliselts.ha=HeatmapAnnotation(CellIslet=cells.tib_2$orig.ident,
                                col=list(CellIslet=islet_colors))

hm6 <- pheatmap(
  af_voi.mat[,cells.tib_2$cell],
  color = color_gradient,
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "AF Heatmap",
  cluster_cols = FALSE,
  top_annotation = celliselts.ha
)

hm7 <- pheatmap(
  af_percentile_mapped[,cells.tib_2$cell],
  color = color_gradient,
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Percentile-Based Heatmap",
  top_annotation = celliselts.ha,
  cluster_cols = FALSE
)
```