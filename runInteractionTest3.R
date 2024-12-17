#!/usr/bin/env Rscript
############################################################
# R script for performing interaction test for pairwise
# 
# Usage: runInteractionTest3.R 
############################################################
knitr::opts_chunk$set(echo = TRUE)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(mashr))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lmerTest))
`%nin%` = Negate(`%in%`)
args = commandArgs(trailingOnly=TRUE) 
############################################################

expr_gene_file <- args[1]
print(expr_gene_file)

geno_snp_file <- args[2]
print(geno_snp_file)

covs_comb_file <- args[3]
print(covs_comb_file)

out_file <- args[4]
print(out_file)

############################################################

cat("\nReading expr_gene")
expr_gene <- read_parquet(expr_gene_file)
print(dim(expr_gene))

# apply INT 
INT <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)), sd=sd(x))
}
# transform expr_gene
expr_gene_t <- data.frame(apply(expr_gene[,-1], 2, INT))
# reset the column names to match the original (matrix changed some dashes to dots introducing mismatches)
colnames(expr_gene_t) <- colnames(expr_gene[,-1])
expr_gene <- cbind(expr_gene[,1], expr_gene_t)
colnames(expr_gene)[1] <- "sample"

cat("\nReading geno_snp")
geno_snp <- read_parquet(geno_snp_file)
print(dim(geno_snp))

cat("\nReading covs_comb")
cov_all <- read_parquet(covs_comb_file)
print(dim(cov_all))
conds <- unique(cov_all$condition)
cov_all$cond <- ifelse(cov_all$condition == conds[1], 1, 2) 
print(table(cov_all$condition, cov_all$cond))
print(head(cov_all))

complete_cols <- sapply(geno_snp, function(col) all(complete.cases(col)))
geno_snp <- geno_snp[, complete_cols]
print(dim(geno_snp))
expr_gene <- expr_gene[, complete_cols]
print(dim(expr_gene))


# 5. Run interaction test 
covs <- data.frame(lapply(cov_all[, !names(cov_all) %in% c("condition", "sample")], as.numeric))
print(head(covs))
print(dim(covs))
# keep only columns with non-missing values
covs <- covs[ , colSums(is.na(covs))==0]
print(dim(covs))

df <- data.frame(gene_id = as.character(),
                 snp_id = as.character(),
                 beta_int = as.numeric(),
                 se_int = as.numeric(),
                 t_int = as.numeric(),
                 p_lmm_int = as.numeric())
for (i in 2:(ncol(expr_gene))){

  tryCatch({
  lmm = lmer(expr_gene[, i] ~ geno_snp[, i] + as.matrix(covs) 
                + geno_snp[, i]:covs[, c("cond")] + (1|cov_all[,"sample"])) # + (1|cov_all[, "condition"])
  lmm_res <- summary(lmm)
  df[nrow(df)+1,] <- list(colnames(expr_gene)[i],  colnames(geno_snp)[i], 
                          lmm_res$coefficients[nrow(lmm_res$coefficients),1],
                          lmm_res$coefficients[nrow(lmm_res$coefficients),2],
                          lmm_res$coefficients[nrow(lmm_res$coefficients),4],
                          lmm_res$coefficients[nrow(lmm_res$coefficients),5])
  if(i%%100==0){
    print(paste(i, "completed", sep=" "))
  }
  }, error = function(e) {
    # Handle the error (you can customize this part)
    cat(paste("Error in model", i, ":", conditionMessage(e), "\n"))
  })
}

write.table(df, out_file)
#arrow::write.parquet(df, sink = out_file, compressed = "uncompressed")
