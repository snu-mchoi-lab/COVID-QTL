###Adopted from https://satijalab.org/seurat/articles/parsebio_sketch_integration
# Install and load required packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
remotes::install_github("bnprks/BPCells/r")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Rsamtools")
install.packages("harmony")

# Load libraries
library(SeuratDisk)
library(Signac)
library(Azimuth)
library(BPCells)
library(Rsamtools)
library(harmony)
library(future)

# Set global options
options(future.globals.maxSize = 8e+09)  # Set memory limit
options(Seurat.object.assay.version = "v5")

# Define file paths
file.dir <- "yourdirectoryhere"
files.set <- c("yourfilename.h5ad")

# Load and process data
data.list <- list()
metadata.list <- list()

for (file in files.set) {
  path <- paste0(file.dir, file)
  
  # Convert and load BPCells matrices
  data <- open_matrix_anndata_hdf5(path)
  output_dir <- paste0(gsub(".h5ad", "", path), "_BP")
  write_matrix_dir(mat = data, dir = output_dir)
  mat <- open_matrix_dir(dir = output_dir)
  
  # Extract metadata
  metadata.list[[file]] <- LoadH5ADobs(path = path)
  data.list[[file]] <- mat
}

# Name layers
names(data.list) <- "COVID"

# Create Seurat object
meta <- LoadH5ADobs(file.path(file.dir, files.set[1]))
COVID <- CreateSeuratObject(counts = data.list, meta.data = meta)

# Normalize and preprocess data
COVID <- NormalizeData(COVID)
COVID <- FindVariableFeatures(COVID)
COVID <- SketchData(COVID, ncells = 200, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(COVID) <- "sketch"

# Perform dimensionality reduction and clustering
COVID <- ScaleData(COVID)
COVID <- RunPCA(COVID)
COVID <- IntegrateLayers(
  object = COVID, 
  method = HarmonyIntegration, 
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = FALSE, 
  group.by = "sample"
)
COVID <- FindNeighbors(COVID, reduction = "harmony", dims = 1:30)
COVID <- FindClusters(COVID, resolution = 2)
COVID <- RunUMAP(COVID, reduction = "harmony", dims = 1:30, return.model = TRUE)

# Visualize clusters
DimPlot(COVID)

# Run Azimuth mapping
COVID <- RunAzimuth(COVID, reference = "pbmcref")

# Project integration onto full dataset
COVID[["sketch"]] <- split(COVID[["sketch"]], f = COVID$sample)
COVID <- ProjectIntegration(
  COVID = COVID, 
  sketched.assay = "sketch", 
  assay = "RNA", 
  reduction = "harmony"
)
COVID <- ProjectData(
  COVID = COVID, 
  sketched.assay = "sketch", 
  assay = "RNA", 
  sketched.reduction = "harmony.full",
  full.reduction = "harmony.full", 
  dims = 1:30, 
  refdata = list(celltype.full = "predicted.celltype.I1")
)

# Visualize full UMAP
DimPlot(COVID, reduction = "umap.full", group.by = "predicted.celltype.I1", alpha = 0.1)

# Perform final clustering and visualization
COVID <- FindNeighbors(COVID, dims = 1:50)
COVID <- FindClusters(COVID, resolution = 2)
COVID <- RunUMAP(COVID, dims = 1:50, return.model = TRUE)
DimPlot(COVID, label = TRUE, label.size = 3, reduction = "umap") + NoLegend()
