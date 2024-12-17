#!/usr/bin/env Rscript

############################################################
# R script for generating list of eQTLs for interaction test
# 
# Usage: genList_int_test_pair.R <conditions> <output_dir> 
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

conditions <- args[1]
print(conditions)

out_folder <- args[2]
print(out_folder) #'/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/interaction_test3/m1_m2/'

############################################################

# separate out the conditions
conds <- strsplit(conditions,"_")[[1]]
print(conds)
cond1 <- conds[1]
cond2 <- conds[2]

## Note: GEX is available for all genes in all conditions

## Keep variants that have genotype info in the two conditions for all cell types

folders <- list.dirs("/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/plink_to_vcf", recursive = F)
fkeep <- grep(paste0(cond1, "|", cond2), folders, value = TRUE)

keep <- NULL
for (f in fkeep){
  
  file <- list.files(f, pattern="*.vcf.gz$", full.names = T) 
  cat("\nReading",file)
  tmp <- fread(file)
  print(dim(tmp)) # 1.14M - 1.18M 
  
  if (is.null(keep)){
    keep <- tmp[,"ID"]
  }
  else{
    keep <- merge(keep, tmp[,"ID"], by="ID")
  }
  cat("\nOutput file: ", dim(keep))
}
keep_fpath <- paste(out_folder, sprintf("vcf_IDs_in_conds_%s_%s.parquet", cond1, cond2), sep='/')
cat("\nWriting to file: ", keep_fpath)
arrow::write_parquet(keep, sink=keep_fpath, compression = "uncompressed")

## Keep ash/mash sig pairs 

if (cond1=="m2" | cond2=="m2"){
  
  # Test example
  # cond1 <- "m1"
  # cond2 <- "m2"
  # out_folder <- "/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/interaction_test3/m1_m2/"
  
  # if either cond1 or cond2 is 'm2', set cond1 to 'm2' and cond2 to the non-m2 condition
  if (cond2=="m2"){
     cond0 = cond1 
     cond1 = "m2"
     cond2 = cond0 
  }
  cat("\nCondition 1: ", cond1)
  cat("\nCondition 2: ", cond2)
  
  # For M2, keep eQTLs from ash (sig) results present in both conditions (keep) and all CTs (m2 mash list) 
  
  flist <-list.files("/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/ash2_results/wo_comorbidity/m2", pattern="*.parquet_sig0.05$", full.names = T)
  cts <- c('B_', 'CD4T_', 'CD8T_', 'DC_', 'Mono_', 'NK_', 'otherT_')
  #keep <- read_parquet(keep_fpath)
  keep$variant_id <- keep$ID
  keep_m2 <- read_parquet("/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/mash/wo_comorbidity/m2/mash_m2_keep.parquet")
  tmp_keep_m2 <- NULL
  for (i in c(1:7)){
    cat("\nReading", flist[i])
    tmp <- read_parquet(flist[i])
    print(dim(tmp))
    # 1. Keep results present in both conditions 
    tmp_keep <- merge(keep[,-c("ID")], tmp[,c(1,2,7:9,18)], by=c("variant_id"))
    cat("\nAccounting for both conditions: ", dim(tmp_keep))
    # 2. Keep results present in all CTs 
    tmp_keep2 <- merge(keep_m2[,c(1,2)], tmp_keep, by=c("phenotype_id","variant_id"))
    cat("\nAccounting for all CTs: ", dim(tmp_keep2))
    # 3. Keep the top eQTL (one per gene) for each cell type for M2
    tmp_keep3 <- tmp_keep2 %>%
      group_by(phenotype_id) %>%
      slice_min(lfsr, with_ties = F)
    cat("\nKeeping top eQTL per CT", dim(tmp_keep3))
    
    # 4. Determine list of eQTLs to use for interaction test for M2
    if (is.null(tmp_keep_m2)){
      tmp_keep_m2 <- tmp_keep3
    }else {
      tmp_keep_m2 <- rbind(tmp_keep_m2, tmp_keep3)
    }
    cat("\nTotal for M2:", dim(tmp_keep_m2))
    tmp_keep_m2 <-unique(tmp_keep_m2[,c("phenotype_id", "variant_id")]) 
    cat("\nTotal unique eQTLs for M2:", dim(tmp_keep_m2))
    
  }
  m2_list_fpath <- paste(out_folder, "ash_sig_results_all_top_m2.parquet", sep='/')
  arrow::write_parquet(tmp_keep_m2, sink=m2_list_fpath, compression = "uncompressed")
  
  # For the non-M2 condition, keep eQTLs in both conditions (keep) and all CTs (mash sig results)
  
  # 1. Read in mash sig results (i.e. present in all CTs)
  fname <- list.files("/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/mash/wo_comorbidity/sig_results", 
                      pattern=sprintf("%s.parquet", cond2), full.names = T)
  cat("\nReading", fname)
  tmp <- read_parquet(fname)
  cat("\nMash sig results for", cond2, ":", dim(tmp))
  # 2. Keep results in both conditions 
  tmp_keep <- merge(keep, tmp, by=c("variant_id"))
  cat("\nAccounting for both conditions: ", dim(tmp_keep))
  # 3. Keep the top eQTL (one per gene) for each cell type for other condition
  cts <- c('B_', 'CD4T_', 'CD8T_', 'DC_', 'Mono_', 'NK_', 'otherT_')
  tmp_keep_cond2 <- NULL
  # Group by gene ID and get rows with lowest value for each lfsr column (cell type)
  for (ct in cts){
    cat("\nCell type:", ct)
    tmp_keep1 <- tmp_keep %>% select(-ID) %>%
      group_by(phenotype_id) %>%
      slice_min(!!sym(paste0(ct, "lfsr")), with_ties = F)
    cat("\t", dim(tmp_keep1))
    
    if (is.null(tmp_keep_cond2)){
      tmp_keep_cond2 <- tmp_keep1
    }
    else {
      tmp_keep_cond2 <- rbind(tmp_keep_cond2, tmp_keep1)
    }
  }
  cat("\nTotal for", cond2, ":", dim(tmp_keep_cond2))
  tmp_keep_cond2 <-unique(tmp_keep_cond2[,c("phenotype_id", "variant_id")]) 
  cat("\nTotal unique eQTLs for", cond2, ":", dim(tmp_keep_cond2))
  cond2_list_fpath <- paste(out_folder, sprintf("mash_sig_results_all_top_%s.parquet", cond2), sep='/')
  arrow::write_parquet(tmp_keep_cond2, sink=cond2_list_fpath, compression = "uncompressed")
  
  # Combine M2 and other condition to generate list of eQTLs for interaction test
  tmp_keep_all <- rbind(tmp_keep_m2, tmp_keep_cond2) 
  tmp_keep_all <- unique(tmp_keep_all)
  all_list_fpath <- paste(out_folder, "sig_top_eQTLs.parquet", sep='/')
  arrow::write_parquet(tmp_keep_all, sink=all_list_fpath, compression = "uncompressed")
  
} else {
  
  # Test case
  # cond1 <- "m1"
  # cond2 <- "m3"
  # out_folder <- "/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/interaction_test3/m1_m3/"
  
  # For both non-M2, keep eQTLs in both conditions (keep) and all CTs (mash sig results)
  flist <- list.files("/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/mash/wo_comorbidity/sig_results", 
                      pattern="*.parquet", full.names = T)
  fkeep <- grep(paste0(cond1, "|", cond2), flist, value = TRUE)
  conds <- c(cond1, cond2)
  cts <- c('B_', 'CD4T_', 'CD8T_', 'DC_', 'Mono_', 'NK_', 'otherT_')
  
  #keep <- read_parquet(keep_fpath)
  keep$variant_id <- keep$ID
  
  tmp_keep_all <- NULL
  for (i in c(1:2)){
    
    # 1. Read in mash sig results (i.e. present in all CTs) and keep results in both conditions
    cat("\nReading", fkeep[i])
    tmp <- read_parquet(fkeep[i])
    print(dim(tmp))
    tmp_keep <- merge(keep, tmp, by=c("variant_id"))
    print(dim(tmp_keep))
    
    # 2. Keep the top eQTL (one per gene) for each cell type for other condition
    tmp_keep2 <- NULL
    # Group by gene ID and get rows with lowest value for each lfsr column (cell type)
    for (ct in cts){
      cat("\nCell type:", ct)
      tmp_keep1 <- tmp_keep %>% select(-ID) %>%
        group_by(phenotype_id) %>%
        slice_min(!!sym(paste0(ct, "lfsr")), with_ties = F)
      cat("\t", dim(tmp_keep1))
      
      if (is.null(tmp_keep2)){
        tmp_keep2 <- tmp_keep1
      }
      else {
        tmp_keep2 <- rbind(tmp_keep2, tmp_keep1)
      }
    }
    cat("\nTotal for", conds[i], ":", dim(tmp_keep2))
    tmp_keep2 <-unique(tmp_keep2[,c("phenotype_id", "variant_id")]) 
    cat("\nTotal unique eQTLs for", conds[i], ":", dim(tmp_keep2))
    cond_list_fpath <- paste(out_folder, sprintf("mash_sig_results_all_top_%s.parquet", conds[i]), sep='/')
    arrow::write_parquet(tmp_keep2, sink=cond_list_fpath, compression = "uncompressed")
    
    if (is.null(tmp_keep_all)){
      tmp_keep_all <- tmp_keep2
    }
    else {
      tmp_keep_all <- rbind(tmp_keep_all, tmp_keep2)
    }
    #print(dim(tmp_keep_all))
  }
  
  # Combine the two conditions to generate list of eQTLs for interaction test
  cat("\nTotal for all conditions:", dim(tmp_keep_all))
  tmp_keep_all <-unique(tmp_keep_all[,c("phenotype_id", "variant_id")]) 
  cat("\nTotal unique eQTLs for all conditions:", dim(tmp_keep_all))
  all_list_fpath <- paste(out_folder, "sig_top_eQTLs.parquet", sep='/')
  arrow::write_parquet(tmp_keep_all, sink=all_list_fpath, compression = "uncompressed")
  
}




