# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("scater", "apeglm"))
remotes::install_github("bnprks/BPCells/r")

# Load libraries
library(BPCells)
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

# Load annotation file
cr38 <- read.table(
  "cellranger_grch38_2020_A_annotations.bed.gz", 
  skip = 5, sep = '\t'
)
colnames(cr38) <- c("#Chr", "start", "end", "ID")
cr38_processed <- cr38

# Define cell types
celltypes <- c("CD4 T", "CD8 T", "B", "DC", "Mono", "NK", "other T")

# Process each cell type
for (celltype in celltypes) {
  CD4_s1 <- subset(COVID_subset, subset = (celltype.full == celltype))
  vcf_samples_order <- read.table(
    "yourfamfile.fam"
  )$V1
  
  counts <- CD4_s1[["RNA"]]$counts
  metadata <- CD4_s1@meta.data
  
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
  
  # Aggregate counts per sample
  groups <- colData(sce)[, "sample"]
  pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "mean")
  tpb <- t(pb)
  
  # Normalize data
  biostars.int <- function(x) { qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))) }
  pseudobulk_int <- apply(tpb, 1, biostars.int) %>% as.data.frame()
  tpseudobulk_int <- pseudobulk_int %>% t() %>% as.data.frame() %>% rownames_to_column(var = "gene")
  
  merged <- merge(cr38_processed, tpseudobulk_int, by.x = "ID", by.y = "gene")
  
  # Reorder and filter columns
  sample_order <- vcf_samples_order[vcf_samples_order %in% colnames(tpseudobulk_int)]
  phenotype_colnames <- c("#Chr", "start", "end", "ID", sample_order)
  
  pheno.data <- merged[, phenotype_colnames]
  pheno.data <- pheno.data[!duplicated(pheno.data$ID), ]
  pheno.data$`#Chr` <- gsub('chr', '', pheno.data$`#Chr`)
  pheno.data <- pheno.data[pheno.data$`#Chr` %in% 1:22, ]
  
  # Save phenotype data
  write.table(
    pheno.data, 
    file = paste0("COVID_", celltype, ".bed"), 
    sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE
  )
}
