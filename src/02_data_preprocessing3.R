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
  rename(ID = Name, MAPINFO = pos) -> HM450K_annotation
HM450_beta %>%
  left_join(HM450K_annotation, by = "ID") %>%
  select(chr, ID, MAPINFO, everything()) %>%
  arrange(chr, MAPINFO) %>%
  filter(chr %in% autosomes) -> HM450_beta

# write beta per chr
lapply(unique(HM450_beta$chr), function(CHR){
  print(CHR)
  HM450_beta[chr == CHR] %>%
    arrange(MAPINFO) %>%
    data.table::fwrite(file = sprintf("tmp/AliciaKSmith_450K_EPIC/HM450_%s.csv", CHR))
})

colnames(EPIC_annotation)
EPIC_annotation %>%
  select(chr, Name, pos) %>%
  rename(ID = Name, MAPINFO = pos) -> EPIC_annotation
EPIC_beta %>%
  left_join(EPIC_annotation, by = "ID") %>%
  select(chr, ID, MAPINFO, everything()) %>%
  arrange(chr, MAPINFO) %>%
  filter(chr %in% autosomes) -> EPIC_beta

# save the data
lapply(unique(EPIC_beta$chr), function(CHR){
  print(CHR)
  EPIC_beta[chr == CHR] %>%
    arrange(MAPINFO) %>%
    data.table::fwrite(file = sprintf("tmp/AliciaKSmith_450K_EPIC/EPIC_%s.csv", CHR))
})


# write training, test, val data of HM450, for the use of PFR
HM450_beta %>%
  select(chr, ID, MAPINFO, all_of(train_samples)) %>%
  arrange(chr, MAPINFO) %>%
  data.table::fwrite("tmp/AliciaKSmith_450K_EPIC/HM450_train.csv")

HM450_beta %>%
  select(chr, ID, MAPINFO, all_of(val_samples)) %>%
  arrange(chr, MAPINFO) %>%
  data.table::fwrite("tmp/AliciaKSmith_450K_EPIC/HM450_val.csv")

HM450_beta %>%
  select(chr, ID, MAPINFO, all_of(test_samples)) %>%
  arrange(chr, MAPINFO) %>%
  data.table::fwrite("tmp/AliciaKSmith_450K_EPIC/HM450_test.csv")







