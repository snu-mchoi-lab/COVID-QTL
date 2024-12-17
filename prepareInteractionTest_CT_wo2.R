#!/usr/bin/env Rscript

############################################################
# R script for performing interaction test
# 
# Usage: prepareInteractionTest_CT_wo2.R <eqtl_list> <ct> <conditions> <out_folder>
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

eqtl_list <- args[1]
print(eqtl_list)

ct <- args[2]
print(ct)

conditions <- args[3]
print(conditions)

out_folder <- args[4]
print(out_folder)

############################################################

# Files to be outputted (combined for conditions)
expr_comb <- NULL
geno_comb <- NULL
covs_comb <- NULL

# Obtain list of conditions 
conds <- strsplit(conditions,"_")[[1]]
print(conds)

# looping through conditions 

for (cond in conds){

    # Read in list of top eQTLs 
    cat("\nReading list of top eQTLs", eqtl_list)
    top <- read_parquet(eqtl_list)

    # Read in gene expression data
    expr_file <- sprintf("/dartfs/rc/lab/S/Szhao/mc_covid_qtl_new/QN_count_matrix/240117_COVID_%s%s_sorted.bed", ct, cond)
    cat("\nReading in gene expression", expr_file)
    expr <- fread(expr_file,header=T)
    print(dim(expr))

    # Read in variant genotype data
    if((ct == "DC_" & cond %in% c("m1", "m3", "s1", "s2", "s3")) | (ct == "otherT_") | (ct == "NK_" & cond %in% c("s3")) ) {
        
        geno_file <- sprintf("/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/plink_to_vcf/%s.%s/covid19_wgs_20230508.qc6.%s.%s.maf0.05.v2.clumpr0.9.final_genotype.vcf.gz", cond, sub("_$", "", ct), cond, sub("_$", "", ct))
        cat("\nReading in variant genotype", geno_file)
        geno <- fread(geno_file, header=T)
        print(dim(geno))

        covariates_file <- sprintf("/dartfs/rc/lab/S/Szhao/mc_covid_qtl_new/Covariates/covid19_wgs_20230508.qc6.%s.%s.maf0.05.v2.clumpr0.9.final.covariate.txt", cond, sub("_$", "", ct)) 

    } else {
        
        geno_file <- sprintf("/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/plink_to_vcf/%s/covid19_wgs_20230508.qc6.%s.maf0.05.v2.clumpr0.9.final_genotype.vcf.gz", cond, cond)
        cat("\nReading in variant genotype", geno_file)
        geno <- fread(geno_file, header=T)
        print(dim(geno))

        covariates_file <- sprintf("/dartfs/rc/lab/S/Szhao/mc_covid_qtl_new/Covariates/230917_covid19_wgs_%s.maf0.05.clump0.9.%s.covariate.txt", cond, sub("_$", "", ct))
    }

    eQTL <- top[, c("phenotype_id", "variant_id")]
    colnames(eQTL) <- c("ID", "variant_id")
    print(dim(eQTL))

    chr_expr <- expr[which(expr$ID %in% c(eQTL$ID)),]
    print(dim(chr_expr))

    chr_geno <- geno[which(geno$ID %in% c(eQTL$variant_id)),] 
    print(dim(chr_geno))

    # Reformat genotype data 
    tmp <- chr_geno[,!c(1:2,4:9)] 
    tmp[tmp=='0/0'] <- 0
    tmp[tmp=='0/1'] <- 1
    tmp[tmp=="1/1"] <- 2
    # Change character genotype to numeric
    tmp1 <- data.frame(tmp)
    cols <- c(names(tmp1)[-1])
    tmp1[cols] <- sapply(tmp1[cols], as.numeric)
    # sapply(tmp1, class)
    chr_snps <- tmp1
    dim(chr_snps)

    # Transform data 
    # gene expression
    tmp <- data.frame(chr_expr[, !c(1:3)])
    #tmp1 <- tmp[,order(names(tmp))]
    expr_trans <- data.frame(t(tmp[,-1]))
    colnames(expr_trans) <- t(tmp[,1])
    expr_trans = tibble::rownames_to_column(expr_trans, "sample")
    print(head(expr_trans)[1:5])
    # snp genotypes
    #tmp <- chr_snps[,order(names(chr_snps))]
    tmp <- chr_snps
    snps_trans <- data.frame(t(tmp[,-1]))
    colnames(snps_trans) <- t(tmp[,1])
    snps_trans = tibble::rownames_to_column(snps_trans, "sample")
    print(head(snps_trans)[1:5])

    # Ensure snp-gene match eQTL pairing
    expr_gene <- expr_trans[, c("sample",eQTL$ID)]
    print(head(expr_gene)[1:6])
    print(dim(expr_gene))
    geno_snp <- snps_trans[, c("sample", eQTL$variant_id)]
    print(head(geno_snp)[1:6])
    print(dim(geno_snp))

    # Add covariates
    cat("\nReading in covariates", covariates_file)
    print(file.exists(covariates_file))
    covariates_comb <- read.table(covariates_file, header=T)
    #covariates_comb <- data.frame(covariates_comb)
    condition <- c("condition", rep(sprintf("%s", cond), ncol(covariates_comb)-1))
    covariates_comb <- rbind(covariates_comb, condition)
    # reformat covariates
    cov_all <- data.frame(t(covariates_comb[,-1]))
    colnames(cov_all) <- t(covariates_comb[,1])
    cov_all$sample = gsub("\\d$", "", rownames(cov_all))
    print(head(cov_all))
    print(dim(cov_all))


    if(is.null(expr_comb)| is.null(geno_comb)| is.null(covs_comb)) {
        expr_comb <- expr_gene 
        geno_comb <- geno_snp
        covs_comb <- cov_all
    }
    else {
        expr_comb <- rbind(expr_comb, expr_gene)
        geno_comb <- rbind(geno_comb, geno_snp)
        covs_comb <- bind_rows(covs_comb, cov_all)

    }

}

cat("\nexpr_comb", dim(expr_comb))
cat("\ngeno_comb", dim(geno_comb))
cat("\ncovs_comb", dim(covs_comb))

arrow::write_parquet(expr_comb, sink=paste(out_folder,sprintf("%sexpr_gene_all_cond.parquet", ct),sep=""))
arrow::write_parquet(geno_comb, sink=paste(out_folder,sprintf("%sgeno_snp_all_cond.parquet", ct),sep=""))
arrow::write_parquet(covs_comb, sink=paste(out_folder,sprintf("%scovs_comb_all_cond.parquet", ct),sep=""))