library(Seurat)
library(tidyverse)
library(Matrix.utils)
library(edgeR)
library(SingleCellExperiment)
library(pheatmap)
library(RColorBrewer)

# Load data
load("yourRData.RData")

# Define cell types
celltypes <- c("CD4_T", "CD8_T", "NK", "B", "Mono", "DC")

# Iterate over each cell type
for (celltype in celltypes) {
  
  # Subset data for the current cell type
  NK <- subset(COVID, subset = (celltype.full == celltype))
  counts <- NK[["RNA"]]$counts
  metadata <- NK@meta.data
  
  # Display metadata summary
  table(metadata$severity, metadata$timepoint)
  
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = metadata
  )
  
  # Aggregate across cluster-sample groups
  groups <- colData(sce)[, "sample"]
  pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "mean")
  tpb <- t(pb)
  
  # Quantile normalization function
  quantile_normalisation <- function(df) {
    df_rank <- apply(df, 2, rank, ties.method = "min")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)
    apply(df_rank, 2, function(idx) df_mean[idx])
  }
  
  # Perform quantile normalization
  pb2 <- quantile_normalisation(tpb)
  write.table(pb2, file = paste0("COVID_", celltype, "_QN.csv"))
  
  # Perform pseudobulk normalization
  biostars.int <- function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
  }
  
  pseudobulk_int <- apply(pb2, 1, biostars.int) %>% as.data.frame()
  tpseudobulk_int <- pseudobulk_int %>% t() %>% as.data.frame() %>% rownames_to_column(var = "gene")
  
  # Merge data
  merged <- merge(cr38_processed, tpseudobulk_int, by.x = "ID", by.y = "gene")
  
  # Define sample order and phenotype column names
  sample_order <- colnames(tpseudobulk_int)
  phenotype_colnames <- c("#Chr", "start", "end", "ID", sample_order)
  
  # Filter phenotype data
  pheno.data <- merged[, phenotype_colnames]
  pheno.data <- pheno.data[!duplicated(pheno.data$ID), ]
  pheno.data$`#Chr` <- gsub('chr', '', pheno.data$`#Chr`)
  pheno.data <- pheno.data[pheno.data$`#Chr` %in% 1:22, ]
  
  # Perform PCA
  expr <- t(pheno.data[, -(1:4)])
  prcompResult <- prcomp(expr, center = TRUE, scale. = FALSE)
  PCs <- prcompResult$x
  write.table(PCs, file = paste0("COVID_", celltype, "_PCs.csv"))
  
  # Determine optimal number of PCs
  resultRunElbow <- PCAForQTL::runElbow(prcompResult = prcompResult)
  print(resultRunElbow)
  write.table(resultRunElbow, file = paste0("COVID_", celltype, "_nPCs.csv"))
  
  # Save phenotype data
  write.table(pheno.data, 
              file = paste0("COVID_", celltype, ".bed"),
              sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
}