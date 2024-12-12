## 0. load required packages
library(SCopeLoomR)
library(SCENIC)
library(dplyr)

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

## 2. extract auc of particular cell type and their information
auc_mtx_{celltype} <- sapply(split(rownames(cellInfo[cellInfo$Celltype_Garnett=='{celltype}',]), cellInfo[cellInfo$Celltype_Garnett=='Mono',]$Sample),
                                                function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
auc_mtx_{celltype}_t <- t(auc_mtx_{celltype})
cellannot <- cellInfo[cellInfo$celltype.full=='{celltype}',c('sample','severity')] 
cellannot <- cellannot %>% distinc()
rownames(cellannot) <- cellannot$sample
                             
## 3. run logistic regression for each cell type
auc_mtx_{celltype}_t <- as.data.frame(auc_mtx_{celltype}_t
data <- merge(auc_mtx_{celltype}_t, cellannot, by='sample') 
ncol(data) #378
p <- c()
est <- c()
or <- c()
ci_up <- c()
ci_low <- c()
for(i in 2:377){
  test <- data[,c(1,378,i)]
  test$severity <- factor(test$severity, levels=c('Mild','Severe'))
  model <- glm(severity ~ test[,3],data=test, family = binomial)
  p <- c(p, coef(summary(model))[2,'Pr(>|z|)'])
  est <- c(est,coef(summary(model))[2,'Estimate'] )
  or <- c(or, exp(coef(model))[2])
  ci_up <- c(ci_up,exp(confint(model))[2,2])
  ci_low <- c(ci_low,exp(confint(model))[2,1])
}
log_df <- data.frame(colnames(data)[2:377],p,est,or,ci_up, ci_low)
colnames(log_df) <- c('TF','p','estimate','oddratio','ci_up','ci_low')
log_df$fdr <- p.adjust(log_df_ot$p,method='fdr')
log_df$bonferroni <- p.adjust(log_df_ot$p,method='bonferroni')
                                      
