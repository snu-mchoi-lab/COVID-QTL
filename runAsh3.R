#!/usr/bin/env Rscript

############################################################
# R script for running ash on QTLs
#
# Usage: run_ash3.R <input_folder> <cell_type> <input_file_type> <input_slope> <input_slope_se> <out_file>
############################################################

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(arrow))
`%nin%` <- Negate(`%in%`)

args <- commandArgs(trailingOnly = TRUE)

############################################################

input_folder <- args[1]
cat("input folder:", input_folder)

cell_type <- args[2]
cat("\ncell type:", cell_type)

input_file_type <- args[3]
cat("\ninput file type:", input_file_type)

input_slope <- args[4]
cat("\nslope:", input_slope)

input_slope_se <- args[5]
cat("\nslope se:", input_slope_se)

out_file <- args[6]
cat("\noutput file:", out_file)

############################################################

list_input <- list.files(input_folder, pattern = sprintf("%s.*\\.%s", cell_type, input_file_type))

allp <- NULL
# loop here is for looping over each chromosome but works for combined data 
for (ls in list_input){

    print(paste("Reading ", input_folder, ls, sep = "/"))
    if (input_file_type == "txt") {
        input1 <- read.table(paste(input_folder, ls, sep = "/"))
    } else if (input_file_type == "parquet"){
        input1 <- read_parquet(paste(input_folder, ls, sep = "/"))
    } else {
        stop("Check file format")
    }

    input1 <- data.frame(input1)
    input1 <- input1[,c(input_slope,input_slope_se)]
    colnames(input1) <- c("slope", "slope_se")
    print(dim(input1))

    if(is.null(allp)){
        allp <- input1
    } else {
        allp <- rbind(allp, input1)
    }
}

cat("\nDimensions of all pairs:", dim(allp))
# Run ash on all pairs 
ash_allp <- ash(allp$slope, allp$slope_se)
print(head(ash_allp$result))

# Add ash results to reference file with all pairs that ash was applied on 
tmp <- cbind(allp, ash_allp$result)
# output ash results file 
arrow::write_parquet(tmp, sink = out_file, compression = "uncompressed")




