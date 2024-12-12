## 0. Required packages
library(coloc)
library(parallel)
library(dplyr)
library(data.table)

## 1. make subset of GWAS data set for lowing computation burden with parallel manner
## 1.1. run this for all 6 conditions respectively (mild1, mild2, mild3, severe1, severer2, and severe3)
disease <- fread('HGIr7c2/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv',header=T)
colnames(disease)[1] = 'CHR'
temp <- disease[disease$all_inv_var_meta_p < 1e-05,]
temp <- temp[!is.na(temp$CHR),]
temp$CHR <- as.numeric(temp$CHR)

subset = function(ct){
  raw.eqtl <- data.frame()
  for(i in seq(1,22){
    filename <- paste0('{path_to_our_eQTL_raw_parquet_file_directory}/230917_',ct,'.mild1.cis_qtl_pairs.',i,'.parquet')
    data <- as.data.frame(read_parquet(filename))
    raw.eqtl <- rbind(raw.eqtl, data)
  }
  raw.eqtl <- cbind(raw.eqtl, str_split_fixed(raw.eqtl$variant_id,':',4)[,c(1,2)])
  colnames(raw.eqtl)[c(10,11)] <- c('chr','pos')
  raw.eqtl$chr <- as.numeric(raw.eqtl$chr)
  raw.eqtl$pos <- as.numeric(raw.eqtl$pos)
  raw.eqtl$id <- paste(raw.eqtl$phenotype_id,raw.eqtl$variant_id,sep=':')
  sig.eqtl.filename <- paste0('{path_to_our_eQTL_eQTL_file_directory}/mild1_',ct,'_filtered_eQTL.txt')
  sig.eqtl <- read.table(sig.eqtl.filename, header=T,sep='\t')
  raw.eqtl$col <- ifelse(raw.eqtl$id %in% sig.eqtl$id, 'red','black')
  print(paste0("making subset eqtl data...[m1-",ct,"]"))
  test.gene <- c()
  for(i in 1:nrow(temp)){
    c = temp[i,]$CHR
    p = temp[i,]$POS
    test.gene <- c(test.gene, unique(raw.eqtl[raw.eqtl$col =='red' & raw.eqtl$chr==c & raw.eqtl$pos >= p -250000 & raw.eqtl$pos <= p +250000,]$phenotype_id))
  }

  subset.eqtl <- raw.eqtl[raw.eqtl$phenotype_id %in% test.gene,]

  filename <- paste0('{path_to_subset_of_gwas_data_directory}/',ct,'.mild1.cis_qtl_pairs.subset_HGIr7.c2.txt')
  print(filename)
  write.table(subset.eqtl, file=filename, col.names=T, row.names=F, quote=F, sep='\t')
}

mclapply(c('B','CD4T','CD8T','DC','NK','Mono'),subset,mc.cores=6)


## 2. load parameter and input files
n_eqtl <- {data frame of sample size with condition information in columns and cell type information in rows}
n_gwas <- {number of sample size in GWAS study}
s <- {proportion of samples in GWAS study that are cases}
colnames(disease) <- c('chr','pos','ref','alt','snp','all_meta_N','all_inv_var_meta_beta','all_inv_var_meta_sebeta','all_inv_var_meta_p','all_inv_var_meta_cases','all_inv_var_meta_controls','all_inv_var_meta_effective','all_inv_var_het_p','lmso_inv_var_beta','lmso_inv_var_se','lmso_inv_var_pval','all_meta_AF','rsid')
disease <- disease[!is.na(disease$chr),]

## 3. run color for each celltype and condition with parallel manner
res <- data.frame()
for (c in c('B','CD4T','CD8T','Mono','NK','DC')){
  for( g in c('m1','m3','s1','s2','s3')){
    input_filename <- paste0('{path_to_subset_of_gwas_data_directory}/',ct,'.',g,'.cis_qtl_pairs.subset_HGIr7.c2.txt')
    input <- read.table(input_filename,header=T)
    i=n[rownames(n)==c,colnames(n)==g]
    
    coloc_peak = function(gene){
      peak.df <- data.frame()
      system(paste("echo 'now processing: ", gene," gene'"))
     
      eqtl.test <- input[input$phenotype_id==gene,]
      chr <-eqtl.test$chr[1]
      gwas.test <- disease[disease$chr == eqtl.test$chr[1] & disease$pos >=min(eqtl.test$pos) & disease$pos <= max(eqtl.test$pos),]
      gwas.test <- gwas.test[!is.na(gwas.test$all_inv_var_meta_beta),]
      gwas.test <- gwas.test[!duplicated(gwas.test$snp),]
      
      if (length(intersect(gwas.test$snp, eqtl.test$variant_id)) >= 20){
        eqtl <- list(eqtl.test$variant_id, eqtl.test$pval_nominal, eqtl.test$slope,(eqtl.test$slope_se)^2,eqtl.test$af,'quant',i)
        names(eqtl) <- c('snp','pval','beta','varbeta','MAF','type','N')
        
        gwas <- list(gwas.test$snp,gwas.test$all_inv_var_meta_p,gwas.test$all_inv_var_meta_beta,gwas.test$all_inv_var_meta_sebeta,'cc',n_gwas,s)
        names(gwas) <- c('snp','pval','beta','varbeta','type','N','s')
        
        result <- coloc.abf(dataset1 = gwas, dataset2=eqtl)
        temp.df <- result$results
        temp.df$testgene <- rep(gene, nrow(temp.df))
        if(nrow(peak.df)==0){
          peak.df<- temp.df
        }else{
          peak.df <- rbind(peak.df, temp.df)
        }
        return(peak.df)  
      }
    }
    mcl_result = mclapply(unique(input$phenotype_id),coloc_peak, mc.cores=10)
    result2 = bind_rows(mcl_result, .id='column_label')
    result2$celltype <- rep(c,nrow(result2))
    result2$group <- rep(g, nrow(result2))
    
    res <- rbind(res, result2)
    print(paste0(c,' celltype, ',g,' group is done'))
  }
}
