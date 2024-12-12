## 0. load required library
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)

## 1. load our ieQTL set
ieqtl <- readRDS('{our_ieqtl_data_frame}')

## 2. make fimo input 
half.window=15

ieqtl$pos <- as.numeric(str_split_fixed(ieqtl$snp,':',4)[,2])
ieqtl$chr <- as.numeric(str_split_fixed(ieqtl$snp,':',4)[,1])
ieqtl$ref <- str_split_fixed(ieqtl$snp,':',4)[,3]
ieqtl$alt <- str_split_fixed(ieqtl$snp,':',4)[,4]
ieqtl.tmp <- ieqtl[,c('snp','pos','chr','ref','alt')] %>% distinct()

meme_list <- list()
for(snp_i in 1:nrow(ieqtl.tmp)){
  snp_id = ieqtl.tmp[snp_i,'snp']
  snp_pos_cur = ieqtl.tmp[snp_i,'pos']
  snp_chrom_cur = ieqtl.tmp[snp_i,'chr']
  ref_allele = ieqtl.tmp[snp_i, 'ref']
  alt_allele = ieqtl.tmp[snp_i, 'alt']
  ref_length = str_length(ref_allele)
  
  left_arm = getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0("chr", snp_chrom_cur), start=(snp_pos_cur - half.window), end=(snp_pos_cur-1))
  right_arm = getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0("chr", snp_chrom_cur), start=(snp_pos_cur + ref_length), end=(snp_pos_cur+ ref_length +half.window))
  
  ref_region = paste0(as.character(left_arm), ref_allele, as.character(right_arm))
  alt_region = paste0(as.character(left_arm), alt_allele, as.character(right_arm))
  meme_tmp = c(paste0(">", snp_id, "_ref"), ref_region,
               paste0(">", snp_id, "_alt"), alt_region)
  
  if(snp_i==1){ meme_list <- meme_tmp}
  else{meme_list <- c(meme_list, meme_tmp)}
}
fileConn <- file('{path_to_input_directory}/ieqtl.fa')
writeLines(meme_list, fileConn)
close(fileConn)

## 3. run fimo in docker
##  docker run -it --rm -v {path_to_working_directory}:/mnt/workspace --name {docker_name} memesuite/memesuite:5.5.1 fimo \
      --oc {path_to_output_directory} \
      --verbosity 1 \
      --thresh 1.0e-4 \
      --max-stored-scores 2000000 \
      /mnt/workspace/{path_to_input_directory}/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
      /mnt/workspace/{path_to_input_directory}/ieqtl.fa

## 4. process output data
fimores <- read.csv('{path_to_output_directory}/fimo.tsv', sep='\t')
fimores = fimores %>% filter(start < 15) %>% filter (stop>15)
fimores = fimores %>% unique()
fimores = fimores %>% tidyr::separate(., col='sequence_name', sep='_',remove=F,into=c('snp','allele'))

fimores_ref = fimores %>% filter(allele == 'ref') %>% dplyr::select(-c("sequence_name", 'allele', 'motif_alt_id')) 
fimores_alt= fimores %>% filter(allele=='alt') %>% dplyr::select(-c("sequence_name", 'allele', 'motif_alt_id')) 
fimores_m = merge(fimores_ref, fimores_alt, by=c("motif_id", 'snp', 'start', 'strand'), suffixes=c("_ref", "_alt"), all=T)
fimores_m = merge(fimores_ref, fimores_alt, by=c("motif_id", 'snp','strand'), suffixes=c("_ref", "_alt"), all=T) 
fimores_m = fimores_m %>% mutate("ref_bind" = (!is.na(`q.value_ref`)) & (`p.value_ref` < 1e-4),
                                 "alt_bind" = (!is.na(`q.value_alt`)) & (`p.value_alt` < 1e-4),
                                 "differential_bind" = ((ref_bind + alt_bind) == 1) )
fimores_d = fimores_m %>% filter(differential_bind)
fimores_d <- fimores_d %>% distinct()
fimores_d$tf <- str_split_fixed(fimores_d$motif_id,'_',2)[,1]


