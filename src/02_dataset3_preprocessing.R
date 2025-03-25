# dataset 3, EPIC V1 and EPIC V2

library(data.table)
library(tidyverse)

V1 <- readRDS("data/EPIC_V1_V2/VHAS_EPICV1_funnorm_combat_n48.rds")
V2 <- readRDS("data/EPIC_V1_V2/VHAS_EPICV2_funnorm_combat_mean_all_probes_n48.rds")
meta <- read.csv("data/EPIC_V1_V2/VHAS_sampleinfo_n48.csv")
meta <- meta[!meta$Technical_replicate,]
V1 <- V1[,meta$Chip_Position_EPICV1]
V2 <- V2[,meta$Chip_Position_EPICV2]

V1_annotation <- fread("data/annotations/GPL21145(EPIC)-48548.txt")
V1_annotation %>%
  mutate(chr = CHR) %>%
  select(chr, ID, MAPINFO) %>%
  na.omit() %>%
  filter(chr %in% as.character(1:22)) %>%
  filter(ID %in% rownames(V1)) %>%
  distinct(ID, .keep_all = T) %>%
  arrange(chr, MAPINFO) -> V1_annotation

V2_annotation <- fread("data/annotations/GPL33022(EPICv2)-233683.txt")
V2_annotation %>%
  mutate(chr = CHR) %>%
  select(chr, ID, MAPINFO) %>%
  na.omit() %>%
  mutate(chr = gsub("chr", "", chr),
         ID = str_split(ID, "_") %>% sapply(`[[`, 1)) %>%
  filter(chr %in% as.character(1:22)) %>%
  filter(ID %in% rownames(V2)) %>%
  distinct(ID, .keep_all = T) %>%
  arrange(chr, MAPINFO) -> V2_annotation

V1 <- V1[V1_annotation$ID,]
V2 <- V2[V2_annotation$ID,]
ggvenn::ggvenn(list(EPIC_V1 = rownames(V1), EPIC_V2 = rownames(V2)),
               fill_color = c("white", "white", "white"),
               auto_scale = T, show_percentage = F)

colnames(V1) <- meta$Sample_Name
colnames(V2) <- meta$Sample_Name

# split into train, validation and test
set.seed(123)
all_samples <- meta$Sample_Name
all_samples %>%
  sample(size = .8 * length(.)) -> train_samples
setdiff(all_samples, train_samples) %>%
  sample(size = .5 * length(.)) -> val_samples
setdiff(all_samples, c(train_samples, val_samples)) -> test_samples

write.table(train_samples, "tmp/EPIC_V1_V2/train_samples.txt",
            row.names = F,
            quote = F,
            col.names = F)

write.table(val_samples, "tmp/EPIC_V1_V2/val_samples.txt",
            row.names = F,
            quote = F,
            col.names = F)

write.table(test_samples, "tmp/EPIC_V1_V2/test_samples.txt",
            row.names = F,
            quote = F,
            col.names = F)

V1 <- cbind(V1_annotation, V1)
V2 <- cbind(V2_annotation, V2)

fwrite(V1, "data/EPIC_V1_V2/EPIC_V1_beta.csv")
fwrite(V2, "data/EPIC_V1_V2/EPIC_V2_beta.csv")

# write beta per chr
lapply(unique(V1$chr), function(x){
  print(x)
  V1[chr == x] %>%
    arrange(MAPINFO) %>%
    data.table::fwrite(file = sprintf("tmp/EPIC_V1_V2/EPIC_V1_chr%s.csv", x))
})

lapply(unique(V2$chr), function(x){
  print(x)
  V2[chr == x] %>%
    arrange(MAPINFO) %>%
    data.table::fwrite(file = sprintf("tmp/EPIC_V1_V2/EPIC_V2_chr%s.csv", x))
})

# write training, test, val data of HM450, for the use of PFR
V1 %>%
  select(chr, ID, MAPINFO, all_of(train_samples)) %>%
  arrange(chr, MAPINFO) %>%
  data.table::fwrite(file = "tmp/EPIC_V1_V2/EPIC_V1_train.csv")

V1 %>%
  select(chr, ID, MAPINFO, all_of(val_samples)) %>%
  arrange(chr, MAPINFO) %>%
  data.table::fwrite(file = "tmp/EPIC_V1_V2/EPIC_V1_val.csv")

V1 %>%
  select(chr, ID, MAPINFO, all_of(test_samples)) %>%
  arrange(chr, MAPINFO) %>%
  data.table::fwrite(file = "tmp/EPIC_V1_V2/EPIC_V1_test.csv")


