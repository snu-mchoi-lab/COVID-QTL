
############################################################
# R script for running mash using significant ash results as strong subset
# 
# Usage: runMash_all2.R <input_tcomb_folder> <input_ash_folder> <sig_level> <nRandom> <nPC> <out_keep_file> <out_data> 
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

input_tcomb_folder <- args[1]
print(input_tcomb_folder)

input_ash_folder <- args[2]
print(input_ash_folder)

sig_level <- args[3]
print(sig_level)

nRandom <- args[4]
print(nRandom)

nPC <- args[5]
print(nPC)

out_keep_file <- args[6]
print(out_keep_file)

out_data <- args[7]
print(out_data)

############################################################

# Function to create matrices B_hat and S_hat (inputs for mash)
create_matrix <- function(df, col, mat){
  
  if (is.null(mat)){
    mat <- df[,col]
  } else{
    mat <- cbind(mat, df[,col])
  }
  
  return(mat)
}

############################################################

# Read in tensorqtl results (chrs combined)
list_tcomb <- list.files(input_tcomb_folder, full.names = T) 
print(list_tcomb)

keep <- NULL 
# Identify pairs in all cell types 
for (f in list_tcomb){
  
  # read all tensorqtl results for all cell types
  print("Reading tensorqtl file")
  tmp <- read_parquet(f)
  cat("new file:",dim(tmp),"\n")
  if (is.null(keep)){
    keep <- tmp
  }
  else{
    keep <- merge(keep[,c(1,2)], tmp[,c(1,2)], by=c("phenotype_id", "variant_id"))  
  }
  cat("keep file:",dim(keep),"\n")
}

arrow::write_parquet(keep, sink = out_keep_file, compression = "uncompressed")

B_hat <- NULL
S_hat <- NULL
# Create mash inputs using all cell types
for (f in list_tcomb){
  
  # read all tensorqtl results for all cell types
  print("Reading tensorqtl file")
  tmp <- read_parquet(f)
  tmp_keep <- merge(tmp, keep, by=c("phenotype_id", "variant_id"))
  print(dim(tmp_keep))
  
  B_hat <- create_matrix(tmp_keep, "slope", B_hat)
  print(dim(B_hat))
  S_hat <- create_matrix(tmp_keep, "slope_se", S_hat)
  print(dim(S_hat))

}


# Label with cell types 
col_names <- c("B", "CD4T","CD8T","DC","Mono","NK","otherT")
colnames(B_hat) <- col_names
colnames(S_hat) <- col_names

# Check for missing values in B_hat or S_hat (remove rows if missing)
if (sum(is.na(B_hat)) != 0 | sum(is.na(S_hat)) != 0) {

  missing_Bhat <- which(is.na(B_hat), arr.ind = TRUE)
  missing_Shat <- which(is.na(S_hat), arr.ind = TRUE)
  missing_indices <- rbind(missing_Bhat, missing_Shat)
  print("Missing indices:")
  print(missing_indices)
  B_hat <- as.matrix(B_hat[-c(unique(missing_indices[,1])),])
  S_hat <- as.matrix(S_hat[-c(unique(missing_indices[,1])),])

} else{
    missing_indices <- NULL
}

# mash dataset for all tests 
data <- mash_set_data(B_hat, S_hat)

### Create strong and random subsets 

## create strong subset 
# use ash results from each cell type 
list_ash2 <- list.files(input_ash_folder, full.names = T)

strong <- NULL
for (i in c(1:length(list_ash2))) {
  
    print("Reading ash2 file")
    tmp <- read_parquet(list_ash2[i])
    print(dim(tmp))
    print(head(tmp))
    
    print("Reading tensorqtl file")
    tensor <- read_parquet(list_tcomb[i])
    tmp <- cbind(tensor[,c(1,2)], tmp)
    print(dim(tmp))
    print(head(tmp))
    
    if (!is.null(missing_indices)){
        tmp <- tmp[-c(unique(missing_indices[,1])),]
    }
    tmp_keep <- merge(tmp, keep, by=c("phenotype_id", "variant_id"))
    print(dim(tmp_keep))
  
    # Get significant ash eQTLs
    ash_sig1 <- tmp_keep[which(tmp_keep$lfsr < as.numeric(sig_level)),]
    print(dim(ash_sig1))
    # Get row numbers of significant ash eQTLs
    index <- as.numeric(rownames(ash_sig1))
  
    if(is.null(strong)){
        strong <- index
    }else{
        strong <- c(strong, index)
    }
}

# get indices for the strong subset (significant in at least one cell type)
strong_subset <- unique(strong)
print(length(strong_subset))

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
data.strong = mash_set_data(B_hat[strong_subset,],S_hat[strong_subset,], V=Vhat)

# obtain covariance matrices 
U.c = cov_canonical(data.random)

U.pca = cov_pca(data.strong, as.numeric(nPC))
U.ed = cov_ed(data.strong, U.pca)

# fit mash model on random tests using cov matrices 
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

# compute posterior summaries 
m_all = mash(data, g=get_fitted_g(m), fixg=TRUE)

print(paste("Number of significant mash eQTLs:",length(get_significant_results(m_all)), sep=" "))
save(m, m_all, file = out_data)

