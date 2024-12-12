## 0. load required packages
library(SCopeLoomR)
library(SCENIC)

## 1. load scenic output
loom <- open_loom('{path_to_output_directory}/{pyscenic.output}.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
cellInfo <- get_cell_annotation(loom)
cellInfo$condition <- paste0(cellInfo$Severity,cellInfo$Timepoint)
meta <- as.data.frame(covid.2e5[[]])
meta$cellbarcode <- rownames(meta)
cellInfo$cellbarcode <- rownames(cellInfo)
rownames(cellInfo) <- cellInfo$cellbarcode


