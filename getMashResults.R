#!/usr/bin/env Rscript

############################################################
# R script for obtaining Mash results
# 
# Usage: getMashResults.R <input_folder> <cond1> <output_folder>
############################################################

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

input_folder <- args[1] # input_folder="/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid/mash_all3/"
print(input_folder)

cond1 <- args[2] # cond1 <- "m1"
print(cond1)

output_folder <- args[3]
print(output_folder)

############################################################
# Get mash results per Condition   

# load mash result output
load(paste(input_folder, sprintf("%s/mash_result_strong_PC5.RData", cond1),sep=""))


if (exists("m_all", envir = .GlobalEnv)) {
  # result for all pairs exist 
  print(dim(get_lfsr(m_all)))
  # get mash significant eQTLs
    #sig_results <- get_significant_results(m_all)
    # get lfsr results
    lfsr_results <- get_lfsr(m_all)
    colnames(lfsr_results) <- c("B_lfsr", "CD4T_lfsr", "CD8T_lfsr", "DC_lfsr", "Mono_lfsr", "NK_lfsr", "otherT_lfsr")
        # colnames(lfsr_results) <- sub("_slope$", "_lfsr", colnames(lfsr_results))
    # get mashr estimated effect sizes 
    mash_beta <- get_pm(m_all)
    colnames(mash_beta) <- c("B_pm", "CD4T_pm", "CD8T_pm", "DC_pm", "Mono_pm", "NK_pm", "otherT_pm")
        # colnames(mash_beta) <- sub("_slope$", "_pm", colnames(mash_beta))

    # identify gene and snp IDs of the mash significant eQTLs 
    tkeep1 = paste(input_folder, sprintf("%s/mash_%s_keep.parquet", cond1,cond1),sep="")
    tmp <- read_parquet(tkeep1)
    tmp_mash_sig <- tmp
    # add lfsr results
    tmp_mash_sig <- cbind(tmp_mash_sig, lfsr_results, mash_beta)
    print(head(tmp_mash_sig))
    print(dim(tmp_mash_sig))

    # write results to parquet file 
    arrow::write_parquet(tmp_mash_sig, sink = paste(output_folder,"mash_all_results_",cond1,".parquet", sep=""), compression = "uncompressed")

} else {
  # results for only strong pairs exist 
  print(dim(get_lfsr(m_strong)))
    # get mash significant eQTLs
    #sig_results <- get_significant_results(m_strong)
    #print(length(sig_results))
    # get lfsr results
    lfsr_results <- get_lfsr(m_strong)
    colnames(lfsr_results) <- c("B_lfsr", "CD4T_lfsr", "CD8T_lfsr", "DC_lfsr", "Mono_lfsr", "NK_lfsr", "otherT_lfsr")
        # colnames(lfsr_results) <- sub("_slope$", "_lfsr", colnames(lfsr_results))
    print(dim(lfsr_results))
    # get mashr estimated effect sizes 
    mash_beta <- get_pm(m_strong)
    colnames(mash_beta) <- c("B_pm", "CD4T_pm", "CD8T_pm", "DC_pm", "Mono_pm", "NK_pm", "otherT_pm")
        # colnames(mash_beta) <- sub("_slope$", "_pm", colnames(mash_beta))
    print(dim(mash_beta))

    # identify gene and snp IDs of the mash significant eQTLs 
    # read in all pairs that are kept
    tkeep1 = paste(input_folder, sprintf("%s/mash_%s_keep.parquet", cond1,cond1),sep="")
    tmp <- read_parquet(tkeep1)
    print(dim(tmp))
    # strong subset
    tmp$min_pval <- apply(tmp[c(seq(from=3, to=21, by=3))], 1, min)
    top <- tmp %>% 
    group_by(phenotype_id) %>%
    # slice_min(min_pval, with_ties = T) # 15560 
    slice_min(min_pval, with_ties = F) # 15539
    print(dim(top))


    tmp_mash_sig <- top
    # add lfsr results
    tmp_mash_sig <- cbind(tmp_mash_sig, lfsr_results, mash_beta)
    print(head(tmp_mash_sig))
    print(dim(tmp_mash_sig))

    # write results to parquet file 
    arrow::write_parquet(tmp_mash_sig, sink = paste(output_folder,"mash_all_results_",cond1,"_strong.parquet", sep=""), compression = "uncompressed")

}

