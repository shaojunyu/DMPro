library(tidyverse)
library(data.table)

# usage: CUE.R <chr> <data_dir> <baseline_dir> [imputation_type]
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("usage: CUE.R <chr> <data_dir> <baseline_dir> [imputation_type]")
}
chr <- as.numeric(args[1])
data_dir <- args[2]
baseline_dir <- args[3]

if (length(args) >= 4) {
  imputation_type <- args[4]
} else if (file.exists(sprintf("%s/EPIC_V1_chr%s.csv", data_dir, chr))) {
  imputation_type <- "EPIC_V1_V2"
} else {
  imputation_type <- "HM450_EPIC"
}

# PFR is an optional member of the ensemble; when its baseline was not run, CUE
# picks among KNN/XGBoost/RF only.
pfr_file <- sprintf("%s/pred_pfr_res_chr%s.rds", baseline_dir, chr)
use_pfr <- file.exists(pfr_file)
if (use_pfr) pfr_res <- readRDS(pfr_file)
glm_res <- readRDS(sprintf("%s/pred_glm_res_chr%s.rds", baseline_dir, chr))
knn_res <- readRDS(sprintf("%s/pred_knn_res_chr%s.rds", baseline_dir, chr))
rf_res <- readRDS(sprintf("%s/pred_rf_res_chr%s.rds", baseline_dir, chr))
xgb_res <- readRDS(sprintf("%s/pred_xgb_res_chr%s.rds", baseline_dir, chr))
cat(sprintf("chr%s  type=%s  use_pfr=%s\n", chr, imputation_type, use_pfr))

train_samples <- read.csv(paste0(data_dir, "/train_samples.txt"), header = F)$V1
test_samples <- read.csv(paste0(data_dir, "/test_samples.txt"), header = F)$V1
val_samples <- read.csv(paste0(data_dir, "/val_samples.txt"), header = F)$V1

# read data; EPIC holds the platform being imputed, HM450 the predictor platform
if (imputation_type == "EPIC_V1_V2") {
  EPIC <- fread(sprintf("%s/EPIC_V2_chr%s.csv", data_dir, chr)) %>% drop_na()
  HM450 <- fread(sprintf("%s/EPIC_V1_chr%s.csv", data_dir, chr)) %>% drop_na()
} else {
  EPIC <- fread(sprintf("%s/EPIC_chr%s.csv", data_dir, chr)) %>% drop_na()
  HM450 <- fread(sprintf("%s/HM450_chr%s.csv", data_dir, chr)) %>% drop_na()
}

setkey(EPIC, ID)
setkey(HM450, ID)

# get the epic-specific probes
setdiff(EPIC$ID, HM450$ID) -> epic_specific
shared_probes <- intersect(EPIC$ID, HM450$ID)

# CUE: per probe, pick the single method with the lowest validation RMSE and use
# its prediction on the test samples.
cue_res <- knn_res %>%
  select(probe, sample, type, y, pred_knn) %>%
  left_join(xgb_res %>% select(probe, sample, pred_xgb), by = c("probe", "sample")) %>%
  left_join(rf_res %>% select(probe, sample, pred_rf), by = c("probe", "sample"))
pred_cols <- c(knn = "pred_knn", xgb = "pred_xgb", rf = "pred_rf")
if (use_pfr) {
  cue_res <- cue_res %>%
    left_join(pfr_res %>% select(probe, sample, pred_PFR), by = c("probe", "sample")) %>%
    replace_na(list(pred_PFR = 0))
  pred_cols <- c(pred_cols, PFR = "pred_PFR")
}

cue_res %>%
  filter(probe %in% epic_specific,
         sample %in% val_samples) %>%
  group_by(probe) %>%
  summarise(across(all_of(unname(pred_cols)), ~ sqrt(mean((y - .x)^2)))) -> val_performance
names(val_performance)[-1] <- names(pred_cols)[match(names(val_performance)[-1], pred_cols)]

min_col <- data.frame(
  probe = val_performance$probe,
  min_col = names(pred_cols)[apply(as.matrix(val_performance[, -1]), 1, which.min)]
)

cue_res %>%
  left_join(min_col, by = "probe") %>%
  mutate(pred_cue = as.numeric(
    as.matrix(pick(all_of(unname(pred_cols))))[cbind(seq_len(n()), match(min_col, names(pred_cols)))]
  )) %>%
  select(probe, pred_cue, y, sample, type) -> cue_res

saveRDS(cue_res, sprintf("%s/pred_cue_res_chr%s.rds", baseline_dir, chr))


