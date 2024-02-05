# This script is used to preprocess the raw data from 450k and EPIC arrays.
# For AnkeHuels_405K_EPIC data
# Pull out the data of autosomes and get the beta values of CpG sites.
# The output is a list of data frames, each of which contains the data of one chromosome.
library(data.table)
library(tidyverse)

# only consider autosomes
chr_names <- paste0("chr", 1:22)

# header of the data frame
HM450_header <- fread("data/AnkeHuels_405K_EPIC/HM450_EPIC/colname_450K.txt") %>% colnames()
EPIC_header <- fread("data/AnkeHuels_405K_EPIC/HM450_EPIC/colname_EPIC.txt") %>% colnames()

# valid samples
HM450_header %>%
  grep("AVG_Beta", ., value = T) %>%
  grep("_R.AVG_Beta", ., invert = T, value = T) %>%
  str_remove(".AVG_Beta") -> HM450_samples

EPIC_header %>%
  grep("AVG_Beta", ., value = T) %>%
  str_replace("_Rep1", "") %>%
  grep("_Rep", ., invert = T, value = T) %>%
  str_remove(".AVG_Beta") -> EPIC_samples

setdiff(HM450_samples, EPIC_samples)
intersect(HM450_samples, EPIC_samples) -> common_samples

# methylation data from 450k array (HumanMethylation450 BeadChip)
HM450_data <- lapply(chr_names, function(chr) {
  cat("Reading 450k data of", chr, "\n")
  fread(paste0("data/AnkeHuels_405K_EPIC/HM450_EPIC/HM450_", chr, ".txt")) %>%
    setnames(., HM450_header) %>%
    mutate(chr = chr)
})
names(HM450_data) <- chr_names

# methylation data from EPIC array
EPIC_data <- lapply(chr_names, function(chr) {
  cat("Reading EPIC data of", chr, "\n")
  fread(paste0("data/AnkeHuels_405K_EPIC/HM450_EPIC/EPIC_", chr, ".txt")) %>%
    setnames(., EPIC_header) %>%
    mutate(chr = chr)
})
names(EPIC_data) <- chr_names

# pull out these columns: [common_samples].AVG_Beta, TargetID, MAPINFO, chr
# sort by MAPINFO
lapply(HM450_data, function(df) {
  df %>%
    select(chr, TargetID, MAPINFO, paste0(common_samples, ".AVG_Beta")) %>%
    arrange(MAPINFO)
}) -> HM450_data

# same for EPIC_data
lapply(EPIC_data, function(df) {
  # rename the columns, remove the _Rep1 in EPIC data
  colnames(df) %>%
    str_replace("_Rep1", "") -> new_colnames
  setnames(df, new = new_colnames)
  df %>%
    select(chr, TargetID, MAPINFO, paste0(common_samples, ".AVG_Beta")) %>%
    arrange(MAPINFO)
}) -> EPIC_data

HM450_merged <- rbindlist(HM450_data) %>%
  setnames(., new = c("chr", "ID", "MAPINFO", common_samples))

HM450_merged %>%
# select(-chr, -MAPINFO) %>%
  drop_na() %>%
  setkey(ID) -> HM450_merged

# write beta per chr
lapply(unique(HM450_merged$chr), function(CHR){
  print(CHR)
  HM450_merged[chr == CHR] %>%
    fwrite(file = sprintf("tmp/processed/HM450_%s.csv", CHR))
})

EPIC_merged <- rbindlist(EPIC_data) %>%
  setnames(., new = c("chr", "ID", "MAPINFO", common_samples))
EPIC_merged %>%
  drop_na() %>%
  setkey(ID) -> EPIC_merged
lapply(unique(EPIC_merged$chr), function(CHR){
  print(CHR)
  EPIC_merged[chr == CHR] %>%
    fwrite(file = sprintf("tmp/processed/EPIC_%s.csv", CHR))
})

# split samples into train, val, test
# train: 80%, val: 10%, test: 10%
set.seed(123)
common_samples %>%
  sample(size = .8 * length(.)) -> train_samples
setdiff(common_samples, train_samples) %>%
  sample(size = .5 * length(.)) -> val_samples
setdiff(common_samples, c(train_samples, val_samples)) -> test_samples

# save the sample names
write.table(train_samples, "tmp/processed/train_samples.txt",
          row.names = F,
          quote = F,
          col.names = F)
write.table(val_samples, "tmp/processed/val_samples.txt",
          row.names = F,
          quote = F,
          col.names = F)
write.table(test_samples, "tmp/processed/test_samples.txt",
          row.names = F,
          quote = F,
          col.names = F)

HM450_merged %>%
  select(chr, ID, MAPINFO, all_of(train_samples)) %>%
  fwrite("tmp/processed/HM450_train.csv")

HM450_merged %>%
  select(chr, ID, MAPINFO, all_of(val_samples)) %>%
  fwrite("tmp/processed/HM450_val.csv")

HM450_merged %>%
  select(chr, ID, MAPINFO, all_of(test_samples)) %>%
  fwrite("tmp/processed/HM450_test.csv")

EPIC_merged %>%
  select(chr, ID, MAPINFO, all_of(train_samples)) %>%
  fwrite("tmp/processed/EPIC_train.csv")

EPIC_merged %>%
  select(chr, ID, MAPINFO, all_of(val_samples)) %>%
  fwrite("tmp/processed/EPIC_val.csv")

EPIC_merged %>%
  select(chr, ID, MAPINFO, all_of(test_samples)) %>%
  fwrite("tmp/processed/EPIC_test.csv")

