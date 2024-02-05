library(GEOquery)
library(tidyverse)
library(data.table)

# read annotation first
anno <- fread("data/annotations/GPL21145(EPIC)-48548.txt")
anno$CHR %>% unique()
anno$MAPINFO

# remove rows whose CHR is NA or MAPINFO is NA
anno <- anno[!is.na(anno$CHR) & !is.na(anno$MAPINFO)]

# order by CHR and MAPINFO, make sure the order of CHR is 1-22, X, Y
anno <- anno[order(anno$CHR, anno$MAPINFO)]
unique(anno$CHR)
fwrite(anno, file = "tmp/processed/EPIC_anno.csv")
anno %>%
  group_by(CHR) %>%
  count() %>%
  write.table("tmp/EPIC_chr_size.txt", row.names = F, col.names = F, quote = F,
              sep = "\t")

#################### GSE197678
# library(minfi)
# ?read.metharray.exp
# RGset <- read.metharray.exp(base = "tmp/GEO/EPIC_data/GSE197678/", force = T,
#                             verbose = T)
# # preprocessing
# Mset <- preprocessIllumina(RGset)
# beta <- getBeta(Mset)
# dim(beta)

beta <- readRDS("beta.rds")
class(beta)

dim(anno)
dim(beta)

rownames(beta)
colnames(beta)
sum(is.na(beta))
# beta[anno$ID,]
anno[anno$ID %in% rownames(beta),] -> anno
anno <- anno[order(anno$CHR, anno$MAPINFO)]
unique(anno$CHR)

beta_dt <- data.table(beta)
dim(beta_dt)
str(beta_dt)
beta_dt$ID <- rownames(beta)
setkey(beta_dt, ID)
beta_dt[anno$ID,] -> beta_dt

setkey(beta_dt, ID)
# write beta per chr
lapply(unique(anno$CHR), function(chr){
  print(chr)
  tmp <- anno$ID[anno$CHR == chr] %>% as.character()
  fwrite(beta_dt[tmp],
         file = sprintf("tmp/processed/GSE197678_chr%s.csv", chr))
})
