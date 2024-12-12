## 0. Load required libraries
library(data.table)
library(stringr)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(Matrix.utils)
library(preprocessCore)
library(tibble)

## 1. refine sample list to test
## 1.0. load HLA allele type table about each individuals 
## 1.1. exclude missed or ambiguous alleles and left first 2 fields.
## 1.2. remain only alleles from HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1, HLA-DQB1, and HLA-DRB1.
## 1.3. exclude samples with too small number genese expressed by more than 3 cells.
## 1.4. save HLA allele type table about remained samples in {allele_table} as dataframe format, and list of remained samples in {sample_list}

## 2. refine gene list to test
## 2.0. exclude genes with expressed in less than 3 cells, with 0 variance, and with most frequent values on over 65% of samples
## 2.1. save expression matrix about each condition with list format in {expr}

gt.all <- readRDS({allele_table})
sample_list <- readRDS({sample_list})
expr <- list()
bad_sample <- list()
sce.filt <- readRDS({seurat_object})
for(t in unique(gt.all$type)){
  
  #make subset of single cell matrix
  celltype <- str_split_fixed(t,':',2)[,2]
  tp <- str_sub(str_split_fixed(t,':',2)[,1],-1)
  samples <- paste0(substr(sample_list[[t]],1,3), '-',substr(sample_list[[t]],4,6),'-', substr(sample_list[[t]],7,9),tp)
  sce.subset <- sce.filt[,sce.filt$celltype.full==celltype & sce.filt$sample %in% samples]
  print(paste0(t,' subset is made'))
  
  #filter out genes with less than 3 non-zero cells
  sce.mat <- data.frame(as.matrix(t(counts(sce.subset))))
  sce.mat$cell <- rownames(sce.mat)
  info <- data.frame(rownames(colData(sce.subset)),colData(sce.subset)$sample)
  colnames(info) <- c('cell','sample')
  sce.mat <- merge(sce.mat,info, by='cell', all.x=T )
  gene_list <- list()
  for(i in unique(info$sample)){
    temp <- sce.mat[sce.mat$sample ==i,][,2:ncol(sce.mat)]
    count_0 <- colSums(temp==0)
    cell_n <- nrow(temp)
    gene <- colnames(temp)[2:ncol(sce.mat)][which(count_0 < cell_n -3 )]
    gene_list[[i]] <- gene
  }
  if(length(gene_list[lengths(gene_list)==0]) >0){
    bad_sample[[t]] <- names(gene_list[lengths(gene_list)==0])
    print(paste0('bad samples were detect in ',t))
    for (n in names(gene_list[lengths(gene_list)==0])){
      gene_list[[n]]<-NULL
    }
  }
  genes.intersect <- Reduce(intersect,gene_list)
  temp <- genes.intersect[ !(genes.intersect %in% all.gene)]
  genes.intersect <- genes.intersect[genes.intersect %in% all.gene]
  temp <- gsub('[.]','-',temp)
  genes.intersect <- c(genes.intersect, temp)
  genes.intersect <- genes.intersect[!is.na(genes.intersect)]
  sce.subset <- sce.subset[genes.intersect,]
  print(paste0(t,' genes filter done #1'))
  
  #pseudobulk
  groups <- colData(sce.subset)[,'sample']
  pb <- aggregate.Matrix(t(counts(sce.subset)),groupings = groups, fun = "mean")
  pb.df <- as.data.frame.matrix(pb)
  pb.qn <- normalize.quantiles(as.matrix(pb.df))
  colnames(pb.qn) <- colnames(pb.df)
  rownames(pb.qn) <- rownames(pb.df)
  t.pb.qn <- t(pb.qn)
  pseudobulk_int <- apply(t.pb.qn, 1, biostars.int) %>% as.data.frame()
  tpseudobulk_int<-t(pseudobulk_int) %>% as.data.frame()%>% rownames_to_column(var="gene")
  print(paste0(t,' pseudo bulk done'))
  
  #filter out genes with 0 variance
  genes.var0 <- tpseudobulk_int$gene[which(apply(tpseudobulk_int[,2:ncol(tpseudobulk_int)],1,var)==0)]
  final.expr <- tpseudobulk_int[! (tpseudobulk_int$gene %in% genes.var0),]
  print(paste0(t,' genes filter done #2'))
  
  #filter out genes with most frequent values on over 65% of samples
  n_sample <- length(sample_list[[t]])
  mf <- apply(final.expr[,2:ncol(final.expr)],1,function(x) names(which.max(table(x))))
  n <- apply(final.expr[,2:ncol(final.expr)],1,function(x) max(table(x)))
  df <- data.frame(final.expr$gene, as.numeric(mf),n)
  genes.high.mf <- df[df$n > n_sample*0.65,]$final.expr.gene
  final.expr <- final.expr[!(final.expr$gene %in% genes.high.mf),]
  print(paste0(t,' genes filter done #3'))
  
  #final.expr$type <- rep(t,nrow(final.expr))
  expr[[t]] <- final.expr
  print(paste0(t, ' is done'))
}

## 3. make covariate data
## 3.0. extract covariate data from eQTL covariate file and save as list format in {cov.all}
             
## 4. make map and ped file for plink input
gene.set <- c('HLA-A', 'HLA-B','HLA-C','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DRB1')
expr <- readRDS({expr})
cov.all <- readRDS({cov.all})
gt.plink <- list()
for(ct in names(sample_list)){
  
  tp = str_sub(str_split_fixed(ct,':',2)[,1],-1)
  sample.set <- sample_list[[ct]]
  allele.info <- allele_table[allele_table$sample %in% sample.set & allele_table$variable %in% gene.set,]
  allele.info$id <- factor(allele.info$id, levels=names(sort(table(allele.info$id),decreasing = F)))
  allele.info$sample <- paste0(substr(allele.info$sample,1,3), '-',substr(allele.info$sample,4,6),'-', substr(allele.info$sample,7,9),tp)
  hla.gene <- unique(allele.info$variable)
  samples <- unique(allele.info$sample)
  gt <- data.frame()
  for(i in unique(allele.info$id)){
    id.temp <- c(i)
    for( s in samples){
      if (s %in% allele.info[allele.info$id ==i,]$sample){
        id.temp <- c(id.temp,1)
      }else{
        id.temp <- c(id.temp,0)
      }
    }
    gt<-rbind(gt,id.temp)
  }
  colnames(gt) <- c('hla_allele',samples)
  neworder <- colnames(expr[[ct]])[2:ncol(expr[[ct]])]
  gt.new <- as.data.frame(gt[,1])
  for(s in neworder){
    gt.new <- cbind(gt.new, gt[,s])
  }
  colnames(gt.new) <- c('hla_allele',neworder)
  samples <- neworder
  gt.plink[[ct]] <- gt.new
  gt <- gt.new
  
  # write map file
  gt$gene <- str_split_fixed(gt$hla_allele,':',3)[,1]
  pos_inf <- gene_anno[gene_anno$gene_id %in% c('HLA-A', 'HLA-B','HLA-C','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DRB1'),]
  map<- data.frame()
  for(i in 1:nrow(gt)){
    if(nrow(map)==0){
      map <- rbind(map, c('CHR','SNP','CM','BP'))
      map <- rbind(map, c(6, gt$hla_allele[i], 0, pos_inf[pos_inf$gene_id==gt$gene[i],]$start))
    }else{
      map <- rbind(map, c(6, gt$hla_allele[i], 0, pos_inf[pos_inf$gene_id==gt$gene[i],]$start))
    }
  }
  map <-map[-1,]
  write.table(map,file=paste0('{directory_for_input}/',str_split_fixed(ct,':',2)[1],'_', str_split_fixed(ct,':',2)[2],'_GT.map'), col.names = F, row.names = F, quote=F, sep='\t')
  
  # write ped file
  gt$ref <- rep('A',nrow(gt))
  gt$alt <- rep('B',nrow(gt))
  ped <- data.frame()
  for(s in samples){
    if(nrow(ped) ==0){
      rsid <- c()
      for( v in gt$hla_allele){
        rsid <- c(rsid, rep(v,2))
      }
      ped <- rbind(ped, c('FID','IID','fatherID','motherID','Sex','phenotype'))
      
      temp <- c(s,s,0,0,cov.all[[ct]][cov.all[[ct]]$id==s,]$sex,-9)
      ped <- rbind(ped, temp)
    }else{
      temp <- c(s,s,0,0,cov.all[[ct]][cov.all[[ct]]$id==s,]$sex,-9)
      ped <- rbind(ped, temp)
    }
  }
  ped <- ped[-1,]
  ped2 <- data.frame()
  ped2 <- rbind(ped2, map$X.SNP.)
  for(i in ped$X.FID.){
    temp <- c()
    for (v in map$X.SNP.){
      if (gt[gt$hla_allele==v,colnames(gt)==i]==1){
        temp <- c(temp,'A B')
      }else{
        temp <- c(temp,'A A')
      }
    }
    ped2 <- rbind(ped2,temp)
  }
  colnames(ped2) <- map$X.SNP.
  ped2 <- ped2[-1,]
  ped <- cbind(ped, ped2)
  write.table(ped,file=paste0('{directory_for_input}/',str_split_fixed(ct,':',2)[1],'_', str_split_fixed(ct,':',2)[2],'_GT.ped'), col.names = F, row.names = F, quote=F, sep='\t')
}

## 5. make bed file about phenotype data 
gene_anno <- fread('cellranger_grch38_2020_A_annotations.bed.gz')
gene_anno <- as.data.frame(gene_anno)
for(ct in names(expr)){
  temp <- merge(gene_anno, expr[[ct]], by.x='gene_id', by.y='gene',all.y=T)
  temp <- temp[,c('#chr','start','end','gene_id', colnames(temp)[5:ncol(temp)])]
  temp <- temp[temp[,1] %in% paste0('chr',seq(1,22)),]
  temp[,1] <- as.numeric(str_split_fixed(temp[,1],'chr',2)[,2])
  temp <- temp[!(duplicated(temp$gene_id)),]
  write.table(temp, file=paste0('{directory_for_input}/',str_split_fixed(ct,':',2)[1],'_', str_split_fixed(ct,':',2)[2],'_GE.bed'),sep='\t',col.names=T,row.names=F, quote=F)
}

## 6. make covariate file 
for(ct in names(cov.all)){  
  temp = cov.all[[ct]]
  temp.t = data.frame(t(temp[,2:ncol(temp)]))
  colnames(temp.t) <- temp$id
  write.table(temp.t, file=paste0('{directory_for_input}/',str_split_fixed(ct,':',2)[1],'_',str_split_fixed(ct,':',2)[2],'_covariate.txt'),col.names=T, row.names=T, sep='\t',quote=F)
}


## 7. convert map and ped file to binary plink file throuth plink
## 8. sort bed file of phenotype data through bedtools
## 9. run tensorqtl
