library(GEOquery)
library(tidyverse)
library(data.table)

# read annotation first
anno <- fread("data/annotations/GPL13534(450K)-11288.txt")
anno$CHR %>% unique()
anno$MAPINFO

# remove rows whose CHR is NA or MAPINFO is NA
anno <- anno[!is.na(anno$CHR) & !is.na(anno$MAPINFO)]

# order by CHR and MAPINFO, make sure the order of CHR is 1-22, X, Y
anno <- anno[order(anno$CHR, anno$MAPINFO)]
fwrite(anno, file = "tmp/processed/450k_anno.csv")

anno %>%
  group_by(CHR) %>%
  count() %>%
  write.table("tmp/450K_chr_size.txt", row.names = F, col.names = F, quote = F,
              sep = "\t")

#################### GSE105018
beta <- fread("tmp/GEO/450k_data/GSE105018_NormalisedData.csv",
              # nrows = 1000,
              check.names = T, header = T
              )
beta %>% rename(ID = V1) -> beta
setkey(beta, ID)
beta[anno$ID,] -> beta

setkey(beta, ID)
# write beta per chr
lapply(unique(anno$CHR), function(chr){
  print(chr)
  tmp <- anno$ID[anno$CHR == chr] %>% as.character()
  fwrite(beta[tmp],
         file = sprintf("tmp/processed/GSE105018_chr%s.csv", chr))
})

# fwrite(beta, file = "tmp/processed/GSE105018.csv")

rm(beta)
gc()
# beta_val <- beta %>%
#   select(!V1) %>%
#   as.data.frame()
# 
# rownames(beta_val) <- beta$V1
# sum(is.na(beta_val))
# 
# cat("aug///")
# beta_val <- beta_val[anno$ID,]
# cat("writing///")

# fwrite(beta, file = "tmp/processed/GSE105018.csv")

# q(save = "no")

#################### GSE55763
beta <- fread("tmp/GEO/450k_data/GSE55763_normalized_betas.txt", 
              # nrows = 100,
              check.names = T)

beta <- beta[ID_REF %in% anno$ID]

# colnames of beta are ID_REF, Sample_1_beta, Sample_1_Pval, Sample_2_beta, Sample_2_Pval, ...
# we need to split them into two data frames, one for beta values, one for p values
beta_val <- beta %>%
  select(!contains("Pval") & !ID_REF) %>%
  as.data.table()
p_val <- beta %>%
  select(contains("Pval")) %>%
  as.data.frame()

ncol(beta_val) == ncol(p_val)
sum(is.na(beta_val)) / nrow(beta_val) / ncol(beta_val)
sum(p_val > 0.05) / nrow(beta_val) / ncol(beta_val)
# hist(beta_val[,1])

# remove probes whose p value > 0.05
beta_val[p_val > 0.05] <- NA
sample_names <- colnames(beta_val)
beta_val$ID <- beta$ID_REF
beta_val[,c("ID", sample_names), with = F] -> beta_val

setkey(beta_val, "ID")
# augment beta_val to the same number of rows as anno, 
# populate NA for probes not in beta_val
beta_val[anno$ID,] -> beta_val

setkey(beta_val, ID)
# write beta per chr
lapply(unique(anno$CHR), function(chr){
  print(chr)
  tmp <- anno$ID[anno$CHR == chr] %>% as.character()
  fwrite(beta_val[tmp],
         file = sprintf("tmp/processed/GSE55763_chr%s.csv", chr))
})

# gc()
# sum(is.na(beta_val)) / nrow(beta_val) / ncol(beta_val)
# fwrite(beta_val, file = "tmp/processed/GSE55763.csv")








