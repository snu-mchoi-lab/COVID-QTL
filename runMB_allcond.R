
############################################################
# R script for running mash-baseline
# 
# Usage: runMB_allcond.R <input_file> <nRandom> <nPC> <out_top_file> <out_data> 
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

input_file <- args[1]
print(input_file)

nRandom <- args[2]
print(nRandom)

nPC <- args[3]
print(nPC)

out_top_file <- args[4]
print(out_top_file)

out_data <- args[5]
print(out_data)

############################################################

set.seed(1234)

# Read in file 
cat("\nReading ", input_file)
keep <- read_parquet(input_file)

# Name the columns with respective conditions
cond <- c("m1", "m3", "s1", "s2", "s3") # m2
for (i in seq(3, ncol(keep), by = 3)) {
  names(keep)[i:(i+2)] <- cond[((i-1)/3 + 1)]
}

# Create Bhat and Shat matrices 
B_hat_cond <- as.matrix(keep[,c(seq(from=4, to=16, by=3))])
print(dim(B_hat_cond))
S_hat_cond <- as.matrix(keep[,c(seq(from=5, to=17, by=3))])
print(dim(S_hat_cond))

# Get indices for a random subset of tests
random_subset <- sample(1:nrow(B_hat_cond), as.numeric(nRandom))
print(length(random_subset)) 

# Obtain correlation structure 
data.temp = mash_set_data(B_hat_cond[random_subset,],S_hat_cond[random_subset,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

# Create data.L to calc z-score
cat("\nCreating mash dataset data_cond")
data_cond = mash_set_data(B_hat_cond,S_hat_cond,V=Vhat)
cat("\nUpdating to create data.L")
data.L = mash_update_data(data_cond, ref = 1)
print(dim(data.L$Bhat))
print(dim(data.L$Shat))
zscore = data.frame(data.L$Bhat/data.L$Shat)
print(dim(zscore))
keepz = cbind(keep[,c(1:2)], zscore)
print(head(keepz))

# Obtain strong subset 
# Top SNPs: largest |z score| across timepoints
get_max <- function(x) {
  abs_row <- abs(x)
  max_abs <- max(abs_row)
  max_value <- x[abs_row == max_abs]
  return(max_value)
}
keepz$max_z <- apply(keepz[,c(3:6)], 1, get_max)
print(head(keepz))
top <- keepz %>% 
  group_by(phenotype_id) %>%
  #slice_max(max_z, with_ties = T) # 15542
  slice_max(max_z, with_ties = F) # 15539
print(dim(top))

# Create B_hat and S_hat matrices for strong subset
keep_top <- merge(top[,c(1,2)], keep, by=c("phenotype_id","variant_id"))
print(dim(keep_top))
arrow::write_parquet(keep_top, sink = out_top_file, compression = "uncompressed")

B_hat_cond_strong <- as.matrix(keep_top[,c(seq(from=4, to=16, by=3))])
print(dim(B_hat_cond_strong))
S_hat_cond_strong <- as.matrix(keep_top[,c(seq(from=5, to=17, by=3))])
print(dim(S_hat_cond_strong))

# Create main data objects for random and strong subsets
data.random = mash_set_data(B_hat_cond[random_subset,],S_hat_cond[random_subset,],V=Vhat)
cat("\nUpdating to create data.random.L")
data.random.L = mash_update_data(data.random, ref = 1)
print(dim(data.random.L$Bhat))
data.strong = mash_set_data(B_hat_cond_strong,S_hat_cond_strong, V=Vhat)
cat("\nUpdating to create data.strong.L")
data.strong.L = mash_update_data(data.strong, ref = 1)
print(dim(data.strong.L$Bhat))

# Obtain covariance matrices 
cat("\nObtaining covariance matrices")
U.c = cov_canonical(data.random.L)
U.pca = cov_pca(data.strong.L, as.numeric(nPC))
U.ed = cov_ed(data.strong.L, U.pca)

# Fitting model using random subset 
cat("\nCalculating model parameters")
m = mash(data.random.L, Ulist = c(U.ed,U.c), outputlevel = 1)
m_pi <- get_estimated_pi(m)
print(m_pi)

# Calculate posterior using model 
cat("\nApplying model to strong subset")
m_strong = mash(data.strong.L, g=get_fitted_g(m), fixg=TRUE)
cat("\nApplying model to all data")
m_all = mash(data.L, g=get_fitted_g(m), fixg=TRUE)

save(m, m_strong, m_all, file = out_data)



