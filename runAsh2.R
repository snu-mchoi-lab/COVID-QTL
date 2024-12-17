#!/usr/bin/env Rscript

############################################################
# R script for running ash on QTLs 
# 
# Usage: run_ash2.R <input_folder> <cell_type> <output_file> 
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

input_folder <- args[1]
print(input_folder)

cell_type <- args[2]
print(cell_type)

out_file <- args[3]
print(out_file)

############################################################

list_tensor <- list.files(input_folder, pattern = cell_type) 

allp <- NULL
# loop here is for looping over each chromosome
for (ls in list_tensor){

    print(paste("Reading ",input_folder,"/ash_",ls,sep=""))
    tensor1 <- read_parquet(paste(input_folder,ls,sep="/"))
    tensor1 <- data.frame(tensor1)
    tensor1 <- tensor1[,c("slope","slope_se")]
    print(dim(tensor1))

    if(is.null(allp)){
        allp <- tensor1
    } else {
        allp <- rbind(allp, tensor1)
    }
}

print(dim(allp))
ash_allp <- ash(allp$slope, allp$slope_se)
print(head(ash_allp$result))

# Add ash results to reference file with all pairs that ash was applied on 
tmp <- cbind(allp, ash_allp$result)
# output ash results file 
arrow::write_parquet(tmp, sink = out_file, compression = "uncompressed")
