## 0. load required packages

## 1. load fimo, scenic, logistic regression, and gene - variant data of each ieQTL
## 1.0. fimo
fimores_d <- readRDS('{fimo_output}')
## 1.1. scenic
loom <- open_loom('{path_to_output_directory}/{pyscenic.output}.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')
embeddings <- get_embeddings(loom)
cellInfo <- get_cell_annotation(loom)
cellInfo$condition <- paste0(cellInfo$Severity,cellInfo$Timepoint)
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Celltype_Garnett),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType.Scaled <- as.data.frame(t(scale(t(regulonActivity_byCellType), center = T, scale=T)))
regulonActivity_byCellType.Scaled$TF <- str_split_fixed(rownames(regulonActivity_byCellType.Scaled),'[(]',2)[,1]
regulonActivity_byCellType.Scaled.melt <- melt(regulonActivity_byCellType.Scaled, id.vars='TF',measure.vars=c('B','CD4T','CD8T','DC','Mono','NK','otherT'))
regulonActivity_byCellType.Scaled.melt$celltype <- regulonActivity_byCellType.Scaled.melt$variable
regulonActivity_byCellType.Scaled.melt$variable <- as.character(regulonActivity_byCellType.Scaled.melt$variable)
regulonActivity_byCellType.Scaled.melt$celltype <- as.character(regulonActivity_byCellType.Scaled.melt$celltype)
## 1.2. merge fimo and scenic 
fimores_d.ct.auc <- merge(fimores_d.ct, regulonActivity_byCellType.Scaled.melt[,c('tf_ct','value')],by='tf_ct',all.x=T)
fimores_d.ct.auc <- fimores_d.ct.auc[!is.na(fimores_d.ct.auc$value),] # value column means scaled AUC
## 1.3. load logistic regresson result and merge it with {fimores_d.ct.auc}
log_df_list <- readRDS('{list of log_df data frame from logistic regression}')
log_df <- data.frame()
for( ct in names(log_df_list)){
  tmp <- log_df_list[[ct]]
  tmp$celltype <- ct
  log_df <- rbind(log_df,tmp)
}
log_df$tf_ct <- paste0(str_split_fixed(log_df$TF,'_',2)[,1],'_',log_df$celltype)
fimores_d.ct.auc.ms <- merge(fimores_d.ct.auc, log_df, by='tf_ct')

## 3. add ieQTL data and ieGene expression data

## 4. add logistic regressoin data

## 5. add eQTL data across the conditions of each ieQTL 
