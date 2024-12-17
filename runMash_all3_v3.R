
############################################################
# R script for obtaining significant ash results and run MASH
# 
# Usage: runMash_all3.R <input_folder> <nRandom> <nPC> <out_keep_file> <out_data> 
############################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(mashr))
`%nin%` = Negate(`%in%`)

args = commandArgs(trailingOnly=TRUE)

############################################################

input_folder <- args[1]
print(input_folder)

nRandom <- args[2]
print(nRandom)

nPC <- args[3]
print(nPC)

out_keep_file <- args[4]
print(out_keep_file)

out_data <- args[5]
print(out_data)

############################################################

# Read in tensorqtl results (chrs combined)
list_tcomb <- list.files(input_folder, pattern="*.txt",full.names = T) 
print(list_tcomb)

keep <- NULL 
# Identify pairs (non-missing values) in all cell types 
for (i in 1:length(list_tcomb)){
  
  # read all results for all cell types
  print("Reading file")
  tmp <- read.table(list_tcomb[i])
  temp <- tmp$gene_id
  tmp$gene_id <- tmp$phenotype_id2
  tmp$phenotype_id2 <- temp
  cat("new file:\n",dim(tmp))
  
  # Check for missing values in slope or slope_se (remove rows if missing)
  missing <- tmp[which(is.na(tmp$beta_int) | is.na(tmp$se_int)),]
  if (nrow(missing) > 0){
    cat("\nmissing values:\n") 
    print(missing)
    tmp <- tmp[which(!is.na(tmp$beta_int) & !is.na(tmp$se_int)),]
    cat("missing removed:\n",dim(tmp))
  }
  
  if (is.null(keep)){
    keep <- tmp[,c(1:4,6)]
  }
  else{
    keep <- merge(keep, tmp[,c(1:4,6)], by=c("gene_id", "snp_id"))  
  }
  cat("\nkeep file:",dim(keep),"\n")
  
}

# rename columns in keep
col_names <- c("B", "CD4T","CD8T","DC","Mono","NK","otherT")
val_names <- c("slope", "slope_se", "p_lmm_int")
colnames(keep) <- c("gene_id", "snp_id",paste(rep(col_names,each=3),rep(val_names,7), sep="_"))

# generate file of eQTL pairs
arrow::write_parquet(keep, sink = out_keep_file, compression = "uncompressed")

#keep <- read_parquet("/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid/interaction_test_all_v2/result/mash_all3_v3/mash_keep.parquet")

# Create B_hat and S_hat matrices 
B_hat <- as.matrix(keep[,c(seq(from=3, to=21, by=3))])
S_hat <- as.matrix(keep[,c(seq(from=4, to=22, by=3))]) 
 
# Create strong subset 
keep$min_pval <- apply(keep[c(seq(from=5, to=23, by=3))], 1, min) 
top <- keep %>% 
  group_by(gene_id) %>%
  # slice_min(min_pval, with_ties = T)  
  slice_min(min_pval, with_ties = F) # 11982
print(dim(top))
# Create B_hat and S_hat matrices 
B_hat_strong <- as.matrix(top[,c(seq(from=3, to=21, by=3))]) 
S_hat_strong <- as.matrix(top[,c(seq(from=4, to=22, by=3))])

# get indices for a random subset of tests
random_subset <- sample(1:nrow(B_hat), as.numeric(nRandom))
print(length(random_subset)) 

### Run mash 

# obtain correlation structure 
data.temp = mash_set_data(B_hat[random_subset,],S_hat[random_subset,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

# create main data objects 
data.random = mash_set_data(B_hat[random_subset,],S_hat[random_subset,],V=Vhat)
data.strong = mash_set_data(B_hat_strong,S_hat_strong, V=Vhat)
data <- mash_set_data(B_hat, S_hat, V=Vhat)

# obtain covariance matrices 
U.c = cov_canonical(data.random)
U.pca = cov_pca(data.strong, as.numeric(nPC))
U.ed = cov_ed(data.strong, U.pca)

# fit mash model on random tests using cov matrices 
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

# compute posterior summaries 
m_strong = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
print(paste("Number of significant mash eQTLs (strong):",length(get_significant_results(m_strong)), sep=" "))

m_all = mash(data, g=get_fitted_g(m), fixg=TRUE)
print(paste("Number of significant mash eQTLs (all):",length(get_significant_results(m_all)), sep=" "))

save(m, m_strong, m_all, file = out_data)

