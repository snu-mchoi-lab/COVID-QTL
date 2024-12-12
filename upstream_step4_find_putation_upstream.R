## 0. load required packages

## 1. load fimo, scenic, logistic regression, and gene - variant data of each ieQTL
fimores_d <- readRDS('{fimo_output}')

regulonActivity_byCellType.Scaled <- as.data.frame(regulonActivity_byCellType.Scaled)
regulonActivity_byCellType.Scaled$TF <- str_split_fixed(rownames(regulonActivity_byCellType.Scaled),'[(]',2)[,1]
regulonActivity_byCellType.Scaled.melt <- melt(regulonActivity_byCellType.Scaled, id.vars='TF',measure.vars=c('B','CD4 T','CD8 T','DC','Mono','NK','other T'))
regulonActivity_byCellType.Scaled.melt$celltype <- regulonActivity_byCellType.Scaled.melt$variable
regulonActivity_byCellType.Scaled.melt$variable <- as.character(regulonActivity_byCellType.Scaled.melt$variable)
regulonActivity_byCellType.Scaled.melt$celltype <- as.character(regulonActivity_byCellType.Scaled.melt$celltype)

## 2. merge with SCENIC data

## 3. add ieQTL data and ieGene expression data

## 4. add logistic regressoin data

## 5. add eQTL data across the conditions of each ieQTL 
