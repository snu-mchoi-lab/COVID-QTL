
############################################################
# R script for obtaining input dataset for runnning mash-baseline 
# 
# Usage: getMB_allcond.R <input_folder> <conditions> <cell_types> <start> <end> 
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

conditions <- args[2]
print(conditions)

cell_types <- args[3]
print(cell_types)

start <- args[4]
print(start)

end <- args[5]
print(end)

############################################################

# Reading all conditions 
# input_folder <- "/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid/mash_all3"
conds <- c("m1", "m3", "s1", "s2", "s3")
#conds <- c(conditions)
print(conds)
ct <- c("B","CD4T","CD8T","DC","Mono","NK","otherT")
#ct <- c(cell_types)
print(ct)
cts <- rep(ct, each=3)
print(cts)

# Get merged datasets 
# start <- 3
# end <- 21
for (i in seq(start, end, by=3)){

    cat("\nCreating merged dataset for cell type: ", cts[i])
    merged <- NULL
    for (c in conds){
        
        cat("\nReading file for condition: ", c)
        tmp <- read_parquet(paste0(input_folder,sprintf("/%s/mash_%s_keep.parquet", c, c)))
        print(head(tmp))
        tmp2 <- tmp[,c(1,2,i:(i+2))]
        print(head(tmp2))

        if (is.null(merged)){
            merged <- tmp2
        }else{
            merged <- merge(merged, tmp2, by=c("phenotype_id","variant_id"))
        }
        print(dim(merged))
    }
    print(dim(merged))
    arrow::write_parquet(merged, sink = sprintf("/dartfs/rc/lab/S/Szhao/liyang/qtl_mapping/covid2/mash_baseline2/wo_comorbidity/merged/%s_allcond_keep.parquet",cts[i]), compression = "uncompressed")
}