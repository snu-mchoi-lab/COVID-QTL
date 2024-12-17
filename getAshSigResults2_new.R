#!/usr/bin/env Rscript

############################################################
# R script for obtaining significant ash results 
# 
# Usage: getAshSigResults2.R <input_ash_folder> <input_tensor_folder> <output_tcomb_folder> <sig_level> <cell_type> <out_folder> 
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

input_ash_folder <- args[1]
print(input_ash_folder)

output_tcomb_folder <- args[2]
print(output_tcomb_folder)

sig_level <- args[3]
print(sig_level)

cell_type <- args[4]
print(cell_type)

out_folder <- args[5]
print(out_folder)


############################################################

list_ash <- list.files(input_ash_folder, pattern = cell_type, full.names = T)

for (a in list_ash){

    # Read ash parquet file 
    print(paste("Reading ",a,sep=""))
    ash_out1 <- read_parquet(a)
    print(summary(ash_out1$lfsr))

    # Get significant ash eQTLs
    ash_sig1 <- ash_out1[which(ash_out1$lfsr < sig_level), ] # !names(ash_out1) %in% c("slope", "slope_se")
    # Get row numbers of significant ash eQTLs
    index <- as.numeric(rownames(ash_sig1))
    print(length(index))
    print(head(index))

    # # Read in tensor parquet files
    # # Combine tensor files for all chromosomes 
    # list_tensor <- list.files(input_tensor_folder, pattern = cell_type, full.names = T) 
    # dfs <- NULL
    # for (t in list_tensor){
    #     df <- read_parquet(t)
    #     if (is.null(dfs)){
    #         dfs <- df
    #     }else{
    #         dfs <- rbind(dfs, df)
    #     }
    # }
    # # Write to output_tcomb_folder
    # arrow::write_parquet(dfs, sink = paste(output_tcomb_folder,"/",basename(a),".parquet", sep=""), compression = "uncompressed")

    # Get info for significant ash eQTLs 
    tcomb <- list.files(output_tcomb_folder, pattern = cell_type, full.names = T) 
    dfs <- read_parquet(tcomb)
    print(dim(dfs))
    tensor_sig1 <- dfs[index,]
    print(dim(tensor_sig1))

    # Combine ash and tensor results for sig eQTLs
    tmp <- cbind(tensor_sig1, ash_sig1)
    print(dim(tmp))

    arrow::write_parquet(tmp, sink = paste(out_folder,"/",basename(a),"_sig",sig_level,sep=''), compression = "uncompressed")

    ## Get top variant per gene 
    # remove duplicate columns of slope and slope_se
    g <- as.integer(ave(names(tmp), names(tmp), FUN = length))
    tmp2 <- tmp[, g == 1 | (rowid(names(tmp)) == 1)]
    top <- tmp2 %>% 
        group_by(phenotype_id) %>% 
        slice_min(lfsr, with_ties = F)
    print(dim(top))

    arrow::write_parquet(top, sink = paste(out_folder,"/",basename(a),"_sig",sig_level,"_top",sep=''), compression = "uncompressed")

}

