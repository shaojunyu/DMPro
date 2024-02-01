library(data.table)
library(eulerr)
library(ggvenn)

# check the overlap between 450k and EPIC
anno_450k <- fread("tmp/processed/450k_anno.csv")
anno_EPIC <- fread("tmp/processed/EPIC_anno.csv")
# 
# # using chr20
# intersect(anno_450k$ID[anno_450k$CHR == "20"],
#           anno_EPIC$ID[anno_EPIC$CHR == "20"])
# 
ggvenn(list("450K_chr20" = anno_450k$ID[anno_450k$CHR == "20"],
            "EPIC_chr20" = anno_EPIC$ID[anno_EPIC$CHR == "20"]),
       text_size = 6,
       fill_color = c("white", "white","white"))

ggvenn(list("450K" = anno_450k$ID,
            "EPIC" = anno_EPIC$ID),
       text_size = 6,
       fill_color = c("white", "white","white"))

tmp <- fread("tmp/processed/450k/GSE55763_chr20.csv")
tmp[,c(2:10)] -> tmp
na.omit(tmp) %>% cor() %>% pheatmap::pheatmap(cluster_rows = F, cluster_cols = F)

na.omit(tmp)
