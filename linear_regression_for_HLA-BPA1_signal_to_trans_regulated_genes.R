## 0. Load required libraries
library(snpStats)
library(data.table)
## 1. load input files
## 1.0. genotype data for the condtion of interest
## 1.1. expression data in DC for the condtion of interest
## 1.2. expression data in CD4T for the condtion of interest

m1.gt <- read.plink('{path_to_plink_directory}/{m1_plink_file}')
m1.gt <- data.frame('GT'=as(m1.gt$genotypes,'numeric')[,c('6:33082722:A:G','6:33011832:T:C')])
m1.gt$sample <- rownames(m1.gt)
colnames(m1.gt) <- c('snp2722','snp1832','sample')
m1.gt$snp2722 <- 2- m1.gt$snp2722
m1.gt$snp1832 <- 2- m1.gt$snp1832

m1.b.exp <- fread('{path_to_expression_matrix_directory}/{m1_B_bed_file}')
m1.b.exp.t <- as.data.frame(t(m1.b.exp))
m1.b.exp.t <- m1.b.exp.t[-c(1,2,3),]
colnames(m1.b.exp.t) <- m1.b.exp.t[1,]
m1.b.exp.t <- m1.b.exp.t[-1,]
m1.b.exp.t$sample <- rownames(m1.b.exp.t)
m1.e.m <- merge(m1.gt, m1.b.exp.t[,c('sample','HLA-DPB1')],by='sample')
m1.e.m$`HLA-DPB1` <- as.numeric(m1.e.m$`HLA-DPB1`)

m1.4t.exp <- fread('{path_to_expression_matrix_directory}/{m1_CD4T_bed_file}')
m1.4t.exp.t <- as.data.frame(t(m1.4t.exp))
m1.4t.exp.t <- m1.4t.exp.t[-c(1,2,3),]
colnames(m1.4t.exp.t) <- m1.4t.exp.t[1,]
m1.4t.exp.t <- m1.4t.exp.t[-1,]
expr <-readRDS('/mnt/workspace/20240314_hla_allele/20240426_expr.rds')
m1.CD4T.filtered.gene <- expr[['mild1:CD4T']]$gene
m1.4t.exp.t <- m1.4t.exp.t[,colnames(m1.4t.exp.t) %in% m1.CD4T.filtered.gene]
m1.4t.exp.t$sample <- rownames(m1.4t.exp.t)
m1.e.m.o <- merge(m1.e.m, m1.4t.exp.t, by='sample')

m1.4t.cov <- fread('20240222_covariates/230917_covid19_wgs_m1.maf0.05.clump0.9.CD4T.covariate.txt')
m1.4t.cov <- as.data.frame(m1.4t.cov)
rownames(m1.4t.cov) <- m1.4t.cov$id
m1.4t.cov <- m1.4t.cov[,-1]
m1.4t.cov.t <- as.data.frame(t(m1.4t.cov))
m1.4t.cov.t$sample <- rownames(m1.4t.cov.t)
m1.e.m.o.c <- merge(m1.e.m.o, m1.4t.cov.t[,c('sample','sex','age','GT_PC1','GT_PC2','GT_PC3')], by='sample')
res.df <- data.frame()
for(i in 5:5717){
  res <- lm(as.numeric(m1.e.m.o.c[,i]) ~ `HLA-DPB1.x`+snp2722 + age+as.factor(sex)+GT_PC1+GT_PC2+GT_PC3,data=m1.e.m.o.c)
  est <- coef(summary(res))[2,1]
  se <- coef(summary(res))[2,2]
  p <- coef(summary(res))[2,4]
  gene <-colnames(m1.e.m.o.c)[i]
  tmp <- data.frame(gene, est, se,p)
  res.df <- rbind(res.df,tmp)
}
res.df$fdr <- p.adjust(res.df$p, method='fdr')
hist(res.df$p)

res.df.snp1832 <- data.frame()
for(i in 5:5717){
  res <- lm(as.numeric(m1.e.m.o.c[,i]) ~ `HLA-DPB1.x`+snp1832 + age+as.factor(sex)+GT_PC1+GT_PC2+GT_PC3,data=m1.e.m.o.c)
  est <- coef(summary(res))[2,1]
  se <- coef(summary(res))[2,2]
  p <- coef(summary(res))[2,4]
  gene <-colnames(m1.e.m.o.c)[i]
  tmp <- data.frame(gene, est, se,p)
  res.df.snp1832 <- rbind(res.df.snp1832,tmp)
}
res.df.snp1832$fdr <- p.adjust(res.df.snp1832$p, method='fdr')
hist(res.df.snp1832$p)#4x7 #20240919_lm_mild1_snp2722_snp1832_pvalue_hist


s1.gt <- read.plink('20240524_rerun_tensorqtl_for_eSNP/covid19_wgs_20230508.qc6.Severe_1.maf0.05.v2.clump0.9.add6_33082722_A_G')
s1.gt <- data.frame('GT'=as(s1.gt$genotypes,'numeric')[,c('6:33082722:A:G','6:33011832:T:C')])
s1.gt$sample <- rownames(s1.gt)
colnames(s1.gt) <- c('snp2722','snp1832','sample')
s1.gt$snp2722 <- 2- s1.gt$snp2722
s1.gt$snp1832 <- 2- s1.gt$snp1832
s1.b.exp <- fread('20240222_gene_expression_matrix/230916_COVID_B_s1.sorted.bed.gz')
s1.b.exp.t <- as.data.frame(t(s1.b.exp))
s1.b.exp.t <- s1.b.exp.t[-c(1,2,3),]
colnames(s1.b.exp.t) <- s1.b.exp.t[1,]
s1.b.exp.t <- s1.b.exp.t[-1,]
s1.b.exp.t$sample <- rownames(s1.b.exp.t)
s1.e.m <- merge(s1.gt, s1.b.exp.t[,c('sample','HLA-DPB1')],by='sample')
s1.e.m$`HLA-DPB1` <- as.numeric(s1.e.m$`HLA-DPB1`)
s1.4t.exp <- fread('20240222_gene_expression_matrix/230916_COVID_CD4T_s1.sorted.bed.gz')
s1.4t.exp.t <- as.data.frame(t(s1.4t.exp))
s1.4t.exp.t <- s1.4t.exp.t[-c(1,2,3),]
colnames(s1.4t.exp.t) <- s1.4t.exp.t[1,]
s1.4t.exp.t <- s1.4t.exp.t[-1,]
expr <-readRDS('/mnt/workspace/20240314_hla_allele/20240426_expr.rds')
s1.CD4T.filtered.gene <- expr[['severe1:CD4T']]$gene
s1.4t.exp.t <- s1.4t.exp.t[,colnames(s1.4t.exp.t) %in% s1.CD4T.filtered.gene]
s1.4t.exp.t$sample <- rownames(s1.4t.exp.t)
s1.e.m.o <- merge(s1.e.m, s1.4t.exp.t, by='sample')

s1.4t.cov <- fread('20240222_covariates/230917_covid19_wgs_s1.maf0.05.clump0.9.CD4T.covariate.txt')
s1.4t.cov <- as.data.frame(s1.4t.cov)
rownames(s1.4t.cov) <- s1.4t.cov$id
s1.4t.cov <- s1.4t.cov[,-1]
s1.4t.cov.t <- as.data.frame(t(s1.4t.cov))
s1.4t.cov.t$sample <- rownames(s1.4t.cov.t)
s1.e.m.o.c <- merge(s1.e.m.o, s1.4t.cov.t[,c('sample','sex','age','GT_PC1','GT_PC2','GT_PC3')], by='sample')

res.df.snp2722.s1 <- data.frame()
for(i in 5:3616){
  res <- lm(as.numeric(s1.e.m.o.c[,i]) ~ `HLA-DPB1`+snp2722 + age+as.factor(sex)+GT_PC1+GT_PC2+GT_PC3,data=s1.e.m.o.c)
  est <- coef(summary(res))[2,1]
  se <- coef(summary(res))[2,2]
  p <- coef(summary(res))[2,4]
  gene <-colnames(s1.e.m.o.c)[i]
  tmp <- data.frame(gene, est, se,p)
  res.df.snp2722.s1 <- rbind(res.df.snp2722.s1,tmp)
}
res.df.snp2722.s1$fdr <- p.adjust(res.df.snp2722.s1$p, method='fdr')
hist(res.df.snp2722.s1$p)

res.df.snp1832.s1 <- data.frame()
for(i in 5:3616){
  res <- lm(as.numeric(s1.e.m.o.c[,i]) ~ `HLA-DPB1`+snp1832 + age+as.factor(sex)+GT_PC1+GT_PC2+GT_PC3,data=s1.e.m.o.c)
  est <- coef(summary(res))[2,1]
  se <- coef(summary(res))[2,2]
  p <- coef(summary(res))[2,4]
  gene <-colnames(s1.e.m.o.c)[i]
  tmp <- data.frame(gene, est, se,p)
  res.df.snp1832.s1 <- rbind(res.df.snp1832.s1,tmp)
}
res.df.snp1832.s1$fdr <- p.adjust(res.df.snp1832.s1$p, method='fdr')
hist(res.df.snp1832.s1$p)#4x7 #20240919_lm_severe1_snp2722_snp1832_pvalue_hist

res.df$type <- 'm1_snp2722'
res.df.snp1832$type <- 'm1_snp1832'
res.df.snp2722.s1$type <- 's1_snp2722'
res.df.snp1832.s1$type <- 's1_snp1832'
res.df <- rbind(res.df, res.df.snp1832)
res.df <- rbind(res.df, res.df.snp2722.s1)
res.df <- rbind(res.df, res.df.snp1832.s1)
