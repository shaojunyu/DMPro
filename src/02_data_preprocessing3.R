library(dplyr)

# read data
HM450_beta <- data.table::fread("data/AliciaKSmith_450K_EPIC/HM450_beta.csv", header = T)
EPIC_beta <- data.table::fread("data/AliciaKSmith_450K_EPIC/EPIC_beta.csv", header = T)
HM450K_annotation <- data.table::fread("data/AliciaKSmith_450K_EPIC/HM450K_annotation.csv")
EPIC_annotation <- data.table::fread("data/AliciaKSmith_450K_EPIC/EPIC_annotation.csv")

colnames(HM450_beta)
colnames(HM450_beta) == colnames(EPIC_beta)

# split samples into train, val, test
# train: 80%, val: 10%, test: 10%
set.seed(123)
common_samples <- intersect(colnames(HM450_beta), colnames(EPIC_beta))
# remove ID
common_samples <- common_samples[common_samples != "ID"]
common_samples
common_samples %>%
  sample(size = .8 * length(.)) -> train_samples
setdiff(common_samples, train_samples) %>%
  sample(size = .5 * length(.)) -> val_samples
setdiff(common_samples, c(train_samples, val_samples)) -> test_samples

write.table(train_samples, "tmp/AliciaKSmith_450K_EPIC/train_samples.txt",
            row.names = F,
            quote = F,
            col.names = F)
write.table(val_samples, "tmp/AliciaKSmith_450K_EPIC/val_samples.txt",
            row.names = F,
            quote = F,
            col.names = F)
write.table(test_samples, "tmp/AliciaKSmith_450K_EPIC/test_samples.txt",
            row.names = F,
            quote = F,
            col.names = F)

colnames(HM450K_annotation)
autosomes <- paste0("chr", 1:22)
HM450K_annotation %>%
  select(chr, Name, pos) %>%
  rename(Chr = chr, ID = Name, MAPINFO = pos) -> HM450K_annotation
HM450_beta %>%
  left_join(HM450K_annotation, by = "ID") %>%
  select(Chr, ID, MAPINFO, everything()) %>%
  arrange(Chr, MAPINFO) %>%
  filter(Chr %in% autosomes) -> HM450_beta

# write beta per chr
lapply(unique(HM450_beta$Chr), function(CHR){
  print(CHR)
  HM450_beta[Chr == CHR] %>%
    arrange(MAPINFO) %>%
    data.table::fwrite(file = sprintf("tmp/AliciaKSmith_450K_EPIC/HM450_%s.csv", CHR))
})

colnames(EPIC_annotation)
EPIC_annotation %>%
  select(chr, Name, pos) %>%
  rename(Chr = chr, ID = Name, MAPINFO = pos) -> EPIC_annotation
EPIC_beta %>%
  left_join(EPIC_annotation, by = "ID") %>%
  select(Chr, ID, MAPINFO, everything()) %>%
  arrange(Chr, MAPINFO) %>%
  filter(Chr %in% autosomes) -> EPIC_beta

# save the data
lapply(unique(EPIC_beta$Chr), function(CHR){
  print(CHR)
  EPIC_beta[Chr == CHR] %>%
    arrange(MAPINFO) %>%
    data.table::fwrite(file = sprintf("tmp/AliciaKSmith_450K_EPIC/EPIC_%s.csv", CHR))
})
