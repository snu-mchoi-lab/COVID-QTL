#!/usr/bin/env Rscript

############################################################
# R script for combining tensorqtl results (combine chromosomes)
#
# Usage: comb_tensorqtl.R <input_tensor_folder> <cell_type> <output_tcomb_file>
############################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(arrow))
`%nin%` = Negate(`%in%`)

args = commandArgs(trailingOnly=TRUE)

############################################################

input_tensor_folder <- args[1]
print(input_tensor_folder)

cell_type <- args[2]
print(cell_type)

output_tcomb_file <- args[3]
print(output_tcomb_file)

############################################################

# Combine tensor files for all chromosomes

# Read in tensor parquet filenames
list_tensor <- list.files(input_tensor_folder, pattern = cell_type, full.names = T)

dfs <- NULL

for (t in list_tensor){

    # read file
    df <- read_parquet(t)

    # combine chrs
    if (is.null(dfs)) {
        dfs <- df
    } else {
        dfs <- rbind(dfs, df)
    }

}

dfs <- data.frame(dfs)
print(dim(dfs))

# Write to output_tcomb_folder
arrow::write_parquet(dfs, sink = output_tcomb_file, compression = "uncompressed")
