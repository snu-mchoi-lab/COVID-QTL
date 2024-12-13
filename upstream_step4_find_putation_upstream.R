## 0. load required packages
library(Seurat)
library(SCopeLoomR)
library(SCENIC)

## 1. load fimo result
fimores_d <- readRDS('{fimo_output}')

## 2. load scenic result
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
                                     
## 3. merge fimo and scenic 
fimores_d.ct.auc <- merge(fimores_d.ct, regulonActivity_byCellType.Scaled.melt[,c('tf_ct','value')],by='tf_ct',all.x=T)
fimores_d.ct.auc <- fimores_d.ct.auc[!is.na(fimores_d.ct.auc$value),] # value column means scaled AUC
                                     
## 4. load logistic regresson result and merge it with {fimores_d.ct.auc}
log_df_list <- readRDS('{list of log_df data frame from logistic regression}')
log_df <- data.frame()
for( ct in names(log_df_list)){
  tmp <- log_df_list[[ct]]
  tmp$celltype <- ct
  log_df <- rbind(log_df,tmp)
}
log_df$tf_ct <- paste0(str_split_fixed(log_df$TF,'_',2)[,1],'_',log_df$celltype)
fimores_d.ct.auc.ms <- merge(fimores_d.ct.auc, log_df, by='tf_ct')

## 5. add ieQTL data and ieGene expression data
## 5.1. load scRNAseq data and generate pseudobulk matrix
## 5.2. run correlation test between each regulon activity and associated ieGene expression
COVID <- readRSD({processed_scRNA_Seurat_object})
COVID.pd <- AggregateExpression(COVID, assays='RNA', return.seurat = TRUE, group.by = c('celltype.full','sample'))

sample.info <- cellInfo[,c('Sample','Severity','condition')]
sample.info <- sample.info %>% distinct()
                                     
cor_i <- c()
for(i in 1:nrow(fimores_d.ct.auc.ms)){
  tf <- fimores_d.ct.auc.ms$tf[i]
  ct <- fimores_d.ct.auc.ms$celltype.x[i]
  iegene<- fimores_d.ct.auc.ms$gene[i]

  # regulon activity
  regulonsAUC.tmp <- regulonsAUC[,rownames(cellInfo[cellInfo$Celltype_Garnett==ct,])]
  cellInfo.tmp <- cellInfo[cellInfo$Celltype_Garnett==ct,]
  regulonActivity_bysample.tmp <- data.frame()
  sample <- c()
  mean <- c()
  for(i in unique(cellInfo.tmp$Sample)){
    sample <- c(sample, i)
    cells <- rownames(cellInfo.tmp[cellInfo.tmp$Sample==i,])
    mean <- c(mean, mean(getAUC(regulonsAUC.tmp[,cells])[rownames(getAUC(regulonsAUC.tmp[,cells]))==paste0(tf,'(+)'),]))
  }
  regulonActivity_bysample.tmp <- mean
  names(regulonActivity_bysample.tmp) <- sample
  tf_sample <- data.frame(names(regulonActivity_bysample.tmp),regulonActivity_bysample.tmp)
  colnames(tf_sample) <- c('sample','tf')
  tf_sample <- merge(tf_sample, sample.info, by.x='sample',by.y='Sample',all.x=T)
  
  # gene expression
  COVID.pd.tmp <- subset(COVID.pd, subset = celltype.full==ct)
  COVID.pd.tmp.data <- GetAssayData(COVID.pd.tmp, slot="data", assay="RNA") 
  COVID.pd.tmp.data <- t(COVID.pd.tmp.data)
  COVID.pd.tmp.data.goi <- as.data.frame(COVID.pd.tmp.data[,iegene])
  rownames(COVID.pd.tmp.data.goi) <- str_split_fixed(rownames(COVID.pd.tmp.data.goi),'_',2)[,2]
  COVID.pd.tmp.data.goi$sample <- rownames(COVID.pd.tmp.data.goi)
  colnames(COVID.pd.tmp.data.goi) <- c('gene','sample')
  tf_sample <- merge(tf_sample, COVID.pd.tmp.data.goi, by='sample',all.x=T)
  cor_i <- c(cor_i,cor(tf_sample$tf, tf_sample$gene))
}
fimores_d.ct.auc.ms$cor <- cor_i

## 5. add eQTL data across the conditions of each ieQTL 
cor_eqtl_beta <- c()
n=0
for(i in 1:nrow(fimores_d.ct.auc.ms)){
  n=n+1
  tf <- fimores_d.ct.auc.ms$tf[i]
  ct <- fimores_d.ct.auc.ms$celltype.x[i]
  iegene<- fimores_d.ct.auc.ms$gene[i]
  iesnp <- fimores_d.ct.auc.ms$snp[i]
  chr <- str_split_fixed(iesnp,':',4)[,1]

  # regulon activity collapse to condition level
  regulonsAUC.tmp <- regulonsAUC[,rownames(cellInfo[cellInfo$Celltype_Garnett==ct,])]
  cellInfo.tmp <- cellInfo[cellInfo$Celltype_Garnett==ct,]
  regulonActivity_bysample.tmp <- data.frame()
  sample <- c()
  mean <- c()
  for(s in unique(cellInfo.tmp$Sample)){
    sample <- c(sample, s)
    cells <- rownames(cellInfo.tmp[cellInfo.tmp$Sample==s,])
    mean <- c(mean, mean(getAUC(regulonsAUC.tmp[,cells])[rownames(getAUC(regulonsAUC.tmp[,cells]))==paste0(tf,'(+)'),]))
  }
  regulonActivity_bysample.tmp <- mean
  names(regulonActivity_bysample.tmp) <- sample
  
  tf_sample <- data.frame(names(regulonActivity_bysample.tmp),regulonActivity_bysample.tmp)
  colnames(tf_sample) <- c('sample','tf')
  tf_sample <- merge(tf_sample, sample.info, by.x='sample',by.y='Sample',all.x=T)
  
  tf_condition_mean <- c(mean(tf_sample[tf_sample$condition=='Mild1',]$tf),mean(tf_sample[tf_sample$condition=='Mild2',]$tf),mean(tf_sample[tf_sample$condition=='Mild3',]$tf),
                         mean(tf_sample[tf_sample$condition=='Severe1',]$tf),mean(tf_sample[tf_sample$condition=='Severe2',]$tf),mean(tf_sample[tf_sample$condition=='Severe3',]$tf))

  # extract  beta value from eqtl data
  m1 <- as.data.frame(read_parquet(paste0('{path_to_tensorqtl_parquet_output}/m1/',ct,'.m1.cis_qtl_pairs.',chr,'.parquet')))
  m2 <- as.data.frame(read_parquet(paste0('{path_to_tensorqtl_parquet_output}/m2/',ct,'.m2.cis_qtl_pairs.',chr,'.parquet')))
  m3 <- as.data.frame(read_parquet(paste0('{path_to_tensorqtl_parquet_output}/m3/',ct,'.m3.cis_qtl_pairs.',chr,'.parquet')))
  s1 <- as.data.frame(read_parquet(paste0('{path_to_tensorqtl_parquet_output}/s1/',ct,'.s1.cis_qtl_pairs.',chr,'.parquet')))
  s2 <- as.data.frame(read_parquet(paste0('{path_to_tensorqtl_parquet_output}/s2/',ct,'.s2.cis_qtl_pairs.',chr,'.parquet')))
  s3 <- as.data.frame(read_parquet(paste0('{path_to_tensorqtl_parquet_output}/s3/',ct,'.s3.cis_qtl_pairs.',chr,'.parquet')))
  
  m1.beta <- m1[m1$phenotype_id==iegene & m1$variant_id==iesnp,]$slope
  if(length(m1.beta)==0){m1.beta <- NA}
  m2.beta <- m2[m2$phenotype_id==iegene & m2$variant_id==iesnp,]$slope
  if(length(m2.beta)==0){m2.beta <- NA}
  m3.beta <- m3[m3$phenotype_id==iegene & m3$variant_id==iesnp,]$slope
  if(length(m3.beta)==0){m3.beta <- NA}
  s1.beta <- s1[s1$phenotype_id==iegene & s1$variant_id==iesnp,]$slope
  if(length(s1.beta)==0){s1.beta <- NA}
  s2.beta <- s2[s2$phenotype_id==iegene & s2$variant_id==iesnp,]$slope
  if(length(s2.beta)==0){s2.beta <- NA}
  s3.beta <- s3[s3$phenotype_id==iegene & s3$variant_id==iesnp,]$slope
  if(length(s3.beta)==0){s3.beta <- NA}
  
  beta <- c(m1.beta,m2.beta, m3.beta,s1.beta,s2.beta, s3.beta)
  
  ## merge
  df_tmp <- data.frame(tf_condition_mean, beta)
  colnames(df_tmp) <- c('regulon_activity_mean','beta')
  df_tmp$severity <- c('mild','mild','mild','severe','severe','severe')
  df_tmp <- df_tmp[!is.na(df_tmp$beta),]
  
  cor_eqtl_beta <- c(cor_eqtl_beta, cor(df_tmp$regulon_activity_mean, abs(df_tmp$beta)))
  print(paste0(n,'th row is done'))
}

## 6. use (1) FIMO result PASS, (2) regulons with abs(log odd raios) >= 50, (3) pairs of regulon and ieGene with abs(Pearson R) >= 0.3 (4) pairs of regulon and eQTL effect size with abs(Pearson R) >= 0.3
