## 0. load required packages
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import anndata as ad
import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
from fbpca import pca
from geosketch import gs

## 1. load our scRNA-seq data
adata= sc.read_h5ad('{our_scRNAdata_raw}.h5ad')

## 2. run geosketch
## 2.0. QC before geosketch
nCells=adata.X.shape[0]
print(nCells)
minCountsPerGene=3*.01*nCells # 3 counts in 1% of cells
print("minCountsPerGene: ", minCountsPerGene)
minSamples=.01*nCells # 1% of cells
print("minSamples: ", minSamples)
sc.pp.filter_cells(adata, min_genes=0)
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
adata = adata[adata.obs['predicted_doublet'] == False, :]
sc.pp.filter_cells(adata, min_genes=200 )
sc.pp.filter_genes(adata, min_cells=3 )
adata = adata[adata.obs['n_genes'] < 4000, :]
adata = adata[adata.obs['percent_mito'] < 0.15, :]
#adata.write('{our_scRNAdata_filtered}.h5ad')

## 2.1. run geosketch
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ["n_counts", "percent_mito"])
sc.pp.scale(adata, max_value=10)
X = adata.X
N = 250000 # Number of samples to obtain from the data set.
sketch_index = gs(X_dimred, N, replace=False)
X_sketch = X_dimred[sketch_index]

## 2.2. extract selective cells from processed scRNA data and save as loom foramt
sketch_cb = adata.obs_names[sketch_index]
adata.scenic= sc.read_h5ad('{our_scRNAdata_raw}.h5ad')
adata.scenic.sketch = adata.scenic[sketch_cb].copy()
f_loom_path_scenic = '{sketch_scRNA_raw}.loom'
row_attrs = {
    "Gene": np.array(adata.scenic.sketch.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.scenic.sketch.obs_names) ,
    "nGene": np.array( np.sum(adata.scenic.sketch.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.scenic.sketch.X.transpose() , axis=0)).flatten() ,
}
lp.create( f_loom_path_scenic, adata.scenic.sketch.X.transpose(), row_attrs, col_attrs)

## 3. run pySENIC in docker
## docker run -it --rm -v {path_to_working_directory}:/mnt/workspace \
##            --name {docker_name} aertslab/pyscenic:0.12.1 pyscenic grn \
##            --num_workers 30 -o /mnt/workspace/{path_to_output_directory}/expr_mat.adjacencies.tsv \
##              /mnt/workspace/{path_to_output_directory}/{sketch_scRNA_raw}.loom \
##              /mnt/workspace/{path_to_provided_reference_by_SCENIC}/allTFs_hg38.txt

## docker run -it --rm \
##    -v {path_to_working_directory}:/mnt/workspace --name {docker_name} aertslab/pyscenic:0.12.1 pyscenic ctx \
##        /mnt/workspace/{path_to_output_directory}/expr_mat.adjacencies.tsv \
##        /mnt/workspace/{path_to_provided_reference_by_SCENIC}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
##        /mnt/workspace/{path_to_provided_reference_by_SCENIC}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
##        --annotations_fname /mnt/workspace/{path_to_provided_reference_by_SCENIC}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
##        --expression_mtx_fname /mnt/workspace/{path_to_output_directory}/{sketch_scRNA_raw}.loom \
##        --mask_dropouts \
##        --output /mnt/workspace/{path_to_output_directory}/{regulons}.csv \
##        --num_workers 35

## docker run -it --rm \
##    -v {path_to_working_directory}:/mnt/workspace --name {docker_name} aertslab/pyscenic:0.12.1 pyscenic aucell \
##        /mnt/workspace/{path_to_output_directory}/{sketch_scRNA_raw}.loom\
##        mnt/workspace/{path_to_output_directory}/{regulons}.csv \
##        -o /mnt/workspace/{path_to_output_directory}/{pyscenic.output}.loom \
##        --num_workers 35

