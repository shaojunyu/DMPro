library(tidyverse)
library(data.table)

# load the training data, chr
# args = commandArgs(trailingOnly=TRUE)
# if (length(args) == 1) {
#   chr <- as.numeric(args[1])
# } else{
#   chr <- 20
# }

pfr_res <- readRDS(sprintf("%s/pred_pfr_res_chr%s.rds", res_dir, chr))
glm_res <- readRDS(sprintf("%s/pred_glm_res_chr%s.rds", res_dir, chr))
knn_res <- readRDS(sprintf("%s/pred_knn_res_chr%s.rds", res_dir, chr))
rf_res <- readRDS(sprintf("%s/pred_rf_res_chr%s.rds", res_dir, chr))
xgb_res <- readRDS(sprintf("%s/pred_xgb_res_chr%s.rds", res_dir, chr))

# train_samples <- read.csv("tmp/processed/train_samples.txt", header = F)$V1
# test_samples <- read.csv("tmp/processed/test_samples.txt", header = F)$V1
# val_samples <- read.csv("tmp/processed/val_samples.txt", header = F)$V1

# read data
# EPIC <- fread(sprintf("tmp/processed/EPIC_chr%s.csv", res_dir, chr)) %>% drop_na()
# HM450 <- fread(sprintf("tmp/processed/HM450_chr%s.csv", res_dir, chr)) %>% drop_na()
# setkey(EPIC, ID)
# setkey(HM450, ID)

# get the epic-specific probes
setdiff(EPIC$ID, HM450$ID) -> epic_specific
shared_probes <- intersect(EPIC$ID, HM450$ID)

# CUE
# 1 / 3 * (knn_res$pred_knn + xgb_res$pred_xgb + rf_res$pred_rf) -> pred_cue
knn_res %>%
  left_join(xgb_res, by = c("probe", "sample")) %>%
  left_join(rf_res, by = c("probe", "sample")) %>%
  left_join(pfr_res, by = c("probe", "sample")) -> cue_res

cue_res %>%
  filter(probe %in% epic_specific,
         sample %in% val_samples) %>%
  select(probe, sample, pred_knn, pred_xgb, pred_rf, pred_PFR, y.x) %>%
  replace_na(list(pred_PFR = 0)) %>%
  group_by(probe) %>%
  summarise(knn = sqrt(mean((y.x - pred_knn)^2)),
            xgb = sqrt(mean((y.x - pred_xgb)^2)),
            rf = sqrt(mean((y.x - pred_rf)^2)),
            PFR = sqrt(mean((y.x - pred_PFR)^2))) -> val_performance
val_performance %>%
  select(-probe) %>%
  # get min value for each row, and the corresponding column name, excluding probe
  mutate(min = apply(., 1, min),
         min_col = apply(., 1, function(x) names(.)[which.min(x)])) %>%
  select(min_col) -> min_col
min_col$probe <- val_performance$probe
cue_res %>%
  left_join(min_col, by = "probe") %>%
  mutate(pred_cue = ifelse(min_col == "knn", pred_knn,
                           ifelse(min_col == "xgb", pred_xgb,
                                  ifelse(min_col == "rf", pred_rf, pred_PFR)))) %>%
  select(probe, pred_cue, y.x, sample, type) %>%
  rename(y = y.x) -> cue_res
saveRDS(cue_res, sprintf("%s/pred_cue_res_chr%s.rds", res_dir, chr))

# cbind(cue_res$pred_knn, cue_res$pred_xgb, cue_res$pred_rf) %>%
#   rowMeans(na.rm = T) -> cue_res$pred_cue


