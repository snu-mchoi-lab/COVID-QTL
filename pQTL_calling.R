# Set working directory
setwd("yourworkingdirectoryhere")

# Load required libraries
library(data.table)
library(dplyr)

# Load expression data
expr_file <- "yourworkingexpressionfilehere"
expr <- fread(expr_file, header = TRUE)
chr_expr <- expr

# Load genotype data
geno_file <- "yourgenotypefilehere.gz"
geno <- fread(geno_file, header = TRUE)
chr_geno <- geno

# Reformat genotype data
tmp <- chr_geno[, -c(1:2, 4:9)]
tmp[tmp == '0/0'] <- 0
tmp[tmp == '0/1'] <- 1
tmp[tmp == '1/1'] <- 2

# Convert genotype data to numeric
tmp1 <- data.frame(tmp)
tmp1[] <- lapply(tmp1, as.numeric)
chr_snps <- tmp1

# Transform expression data
expr_trans <- t(as.matrix(chr_expr[, -c(1:3)]))
colnames(expr_trans) <- chr_expr[[1]]

# Transform SNP genotype data
snps_trans <- t(as.matrix(chr_snps))
colnames(snps_trans) <- chr_snps[[1]]

# Load covariate data
covariates_file <- "yourcovariatefilehere.txt"
covariates_comb <- read.table(covariates_file, header = TRUE)

# Reformat covariates
cov_all <- t(as.matrix(covariates_comb[, -1]))
colnames(cov_all) <- covariates_comb[[1]]
rownames(cov_all) <- gsub("\\.", "-", rownames(cov_all))
cov_all[cov_all == '0/0'] <- 0
cov_all[cov_all == '0/1'] <- 1
cov_all[cov_all == '1/1'] <- 2

# Convert covariate columns to numeric
df <- as.data.frame(cov_all, stringsAsFactors = FALSE)
df <- df %>% mutate(across(everything(), as.numeric))
cov_all_clean <- as.matrix(df)
cov_all <- cov_all_clean

# Function to perform linear regression
perform_regression <- function(snp_column, gene_column, covariates) {
  model <- lm(gene_column ~ snp_column + covariates)
  p_value <- summary(model)$coef[2, 4]  # Extract p-value
  
  if (p_value < 0.05) {
    estimate <- summary(model)$coef[2, 1]
    se <- summary(model)$coef[2, 2]
    t_value <- summary(model)$coef[2, 3]
    return(c(estimate, se, t_value, p_value))
  } else {
    return(NULL)
  }
}

# Open a CSV file for writing
output_file <- "yourresultfilename.csv"
cat("snp_col,gene_col,estimate,se,t_value,p_value\n", file = output_file)

# Iterate through all SNP and gene pairs
for (snp_col in colnames(snps_trans)) {
  for (gene_col in colnames(expr_trans)) {
    cat("Linear Regression for SNP:", snp_col, "and Gene:", gene_col, "\n")
    result <- perform_regression(snps_trans[, snp_col], expr_trans[, gene_col], cov_all)
    
    if (!is.null(result)) {
      cat(snp_col, gene_col, paste(result, collapse = ","), sep = ",", file = output_file, append = TRUE)
      cat("\n", file = output_file, append = TRUE)
    }
  }
}

cat("Analysis complete. Results saved to:", output_file, "\n")