# This script is used to preprocess the raw data from 450k and EPIC arrays.
# Pull out the data of autosomes and get the beta values of CpG sites.
# The output is a list of data frames, each of which contains the data of one chromosome.

library(data.table)
library(tidyverse)

# only consider autosomes
chr_names <- paste0("chr", 1:22)

# header of the data frame
HM450_header <- fread("data/HM450_EPIC/colname_450K.txt") %>% colnames()
EPIC_header <- fread("data/HM450_EPIC/colname_EPIC.txt") %>% colnames()

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
  fread(paste0("data/HM450_EPIC/HM450_", chr, ".txt")) %>%
    setnames(., HM450_header) %>%
    mutate(chr = chr)
})
names(HM450_data) <- chr_names

# methylation data from EPIC array
EPIC_data <- lapply(chr_names, function(chr) {
  cat("Reading EPIC data of", chr, "\n")
  fread(paste0("data/HM450_EPIC/EPIC_", chr, ".txt")) %>%
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
  setnames(., new = c("chr", "TargetID", "MAPINFO", common_samples))
fwrite(HM450_merged, "data/HM450_beta_value_merged.csv")

EPIC_merged <- rbindlist(EPIC_data) %>%
  setnames(., new = c("chr", "TargetID", "MAPINFO", common_samples))
fwrite(EPIC_merged, "data/EPIC_beta_value_merged.csv")


# split samples into train, val, test
# train: 70%, val: 10%, test: 20%
set.seed(123)
common_samples %>%
  sample(size = .7 * length(.)) -> train_samples
setdiff(common_samples, train_samples) %>%
  sample(size = .33 * length(.)) -> val_samples
setdiff(common_samples, c(train_samples, val_samples)) -> test_samples

# save the sample names
write.table(train_samples, "data/train_samples.txt",
          row.names = F,
          quote = F,
          col.names = F)
write.table(val_samples, "data/val_samples.txt",
          row.names = F,
          quote = F,
          col.names = F)
write.table(test_samples, "data/test_samples.txt",
          row.names = F,
          quote = F,
          col.names = F)

HM450_merged %>%
  select(chr, TargetID, MAPINFO, all_of(train_samples)) %>%
  fwrite("tmp/HM450_train.csv")

HM450_merged %>%
  select(chr, TargetID, MAPINFO, all_of(val_samples)) %>%
  fwrite("tmp/HM450_val.csv")

HM450_merged %>%
  select(chr, TargetID, MAPINFO, all_of(test_samples)) %>%
  fwrite("tmp/HM450_test.csv")

EPIC_merged %>%
  select(chr, TargetID, MAPINFO, all_of(train_samples)) %>%
  fwrite("tmp/EPIC_train.csv")

EPIC_merged %>%
  select(chr, TargetID, MAPINFO, all_of(val_samples)) %>%
  fwrite("tmp/EPIC_val.csv")

EPIC_merged %>%
  select(chr, TargetID, MAPINFO, all_of(test_samples)) %>%
  fwrite("tmp/EPIC_test.csv")

