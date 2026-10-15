# Cross-dataset (transfer) baselines for CUE: train RF / XGBoost / KNN / GLM
# (and PFR with RUN_PFR=1) on the SOURCE dataset, predict the TARGET dataset's
# test samples. Mirrors src/CUE/baseline.R.
# usage: Transfer_baseline.R <chr> <source_dir> <target_dir> <res_dir>
library(data.table)
library(tidyverse)
options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
options(ignore.interactive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: Transfer_baseline.R <chr> <source_dir> <target_dir> <res_dir>")
}
chr        <- as.numeric(args[1])
source_dir <- args[2]           # models trained on this dataset's train samples
target_dir <- args[3]           # predictions made for this dataset's test samples
res_dir    <- args[4]

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("chr%s  source=%s  target=%s  res=%s\n",
            chr, source_dir, target_dir, res_dir))

# ---- load & align source / target on shared probes ------------------------
HM450        <- fread(sprintf("%s/HM450_chr%s.csv", source_dir, chr))
EPIC         <- fread(sprintf("%s/EPIC_chr%s.csv",  source_dir, chr))
target_HM450 <- fread(sprintf("%s/HM450_chr%s.csv", target_dir, chr))
target_EPIC  <- fread(sprintf("%s/EPIC_chr%s.csv",  target_dir, chr))

# keep the target's full 450K probe list to define its EPIC-specific probes below
target_hm450_full_ids <- target_HM450$ID

HM450_probes <- intersect(HM450$ID, target_HM450$ID)
EPIC_probes  <- intersect(EPIC$ID,  target_EPIC$ID)

HM450        <- HM450[match(HM450_probes, HM450$ID), ]
EPIC         <- EPIC[match(EPIC_probes, EPIC$ID), ]
target_HM450 <- target_HM450[match(HM450_probes, target_HM450$ID), ]
target_EPIC  <- target_EPIC[match(EPIC_probes, target_EPIC$ID), ]

stopifnot(all(EPIC$ID == target_EPIC$ID), all(HM450$ID == target_HM450$ID))

# guard against sample-name collisions between the two datasets
src_epic_samples <- names(EPIC)[4:ncol(EPIC)]
tgt_epic_samples <- names(target_EPIC)[4:ncol(target_EPIC)]
if (length(intersect(src_epic_samples, tgt_epic_samples)) > 0) {
  stop("source and target share sample names; cbind would collide")
}

# combine sample columns: source samples (for training/val) + target samples (test)
EPIC  <- cbind(EPIC,  target_EPIC[,  4:ncol(target_EPIC)])
HM450 <- cbind(HM450, target_HM450[, 4:ncol(target_HM450)])

HM450 <- HM450 %>% drop_na() %>% distinct(ID, .keep_all = TRUE) %>% arrange(MAPINFO)
EPIC  <- EPIC  %>% drop_na() %>% distinct(ID, .keep_all = TRUE) %>% arrange(MAPINFO)

# imputation target: probes on the target's EPIC but not on its 450K
target_specific <- setdiff(target_EPIC$ID, target_hm450_full_ids)
EPIC <- EPIC %>% filter(ID %in% target_specific)

# ---- sample splits: train/val from SOURCE, test from TARGET ----------------
train_samples <- read.csv(sprintf("%s/train_samples.txt", source_dir), header = F)$V1
val_samples   <- read.csv(sprintf("%s/val_samples.txt",   source_dir), header = F)$V1
test_samples  <- read.csv(sprintf("%s/test_samples.txt",  target_dir), header = F)$V1
train_samples <- intersect(train_samples, names(EPIC))
val_samples   <- intersect(val_samples,   names(EPIC))
test_samples  <- intersect(test_samples,  names(EPIC))
cat(sprintf("train(src)=%d  val(src)=%d  test(tgt)=%d  probes=%d\n",
            length(train_samples), length(val_samples), length(test_samples), nrow(EPIC)))

# ---- for each EPIC probe, gather 25 up/down HM450 neighbours ---------------
EPIC  <- EPIC  %>% arrange(MAPINFO)
HM450 <- HM450 %>% arrange(MAPINFO)

y_probes <- EPIC[, .(ID, MAPINFO)]
predictor_names <- paste0("X", 1:50)
no_neighbours <- setNames(as.data.frame(matrix(rep("X", 50), nrow = 1)), predictor_names)
pbmcapply::pbmclapply(1:nrow(y_probes), function(i){
  y_probes[i, .(ID, MAPINFO)] -> tmp
  dist <- abs(tmp$MAPINFO - HM450$MAPINFO)
  which.min(dist) -> idx
  if (idx <= 24 || idx > (nrow(HM450) - 25)) return(no_neighbours)
  HM450[(idx - 24):(idx + 25), .(ID, MAPINFO)] -> tmp2
  setNames(as.data.frame(matrix(tmp2$ID, nrow = 1)), predictor_names)
}) -> near_probes

rbindlist(near_probes) -> X_probes
stopifnot(identical(names(X_probes), predictor_names))
predictors <- cbind(y_probes, X_probes)
predictors[predictors$X1 != "X"] -> predictors
setkey(predictors, ID)

EPIC <- EPIC %>% filter(ID %in% predictors$ID)
setkey(EPIC, ID)
EPIC <- EPIC[predictors$ID]
stopifnot(all(EPIC$ID == predictors$ID))

setkey(HM450, ID)
generate_data <- function(probe_name, EPIC, HM450, predictors){
  EPIC %>%
    filter(ID == probe_name) %>%
    select(-MAPINFO, -chr) %>%
    pivot_longer(!ID, names_to = "sample_name", values_to = "beta") %>%
    as.data.frame() -> y
  rownames(y) <- y$sample_name

  predictors %>%
    filter(ID == probe_name) %>%
    select(starts_with("X")) %>% unlist() %>% as.character() -> tmp_predictors

  HM450[tmp_predictors, ] %>%
    select(-chr, -MAPINFO) %>%
    select(ID, all_of(y$sample_name)) %>%
    as.data.frame() -> X
  rownames(X) <- X$ID
  X[, -1] -> X
  X <- t(X) %>% as.matrix()
  return(list(X, y))
}

############################# Random Forest ####################################
library(randomForest)
cat("Random Forest\n")
pbmcapply::pbmclapply(EPIC$ID, function(x){
  generate_data(x, EPIC, HM450, predictors) -> tmp
  X <- tmp[[1]]; y <- tmp[[2]]
  randomForest(x = X[train_samples, ], y = y[train_samples, ]$beta,
               ntree = 2, mtry = 3, verbose = 0) -> rf
  pred_RF_test <- predict(rf, X[test_samples, ])
  pred_RF_val  <- predict(rf, X[val_samples, ])
  test_res <- data.frame(probe = x, pred_rf = pred_RF_test, y = y[test_samples, ]$beta,
                         sample = test_samples, type = "test")
  val_res  <- data.frame(probe = x, pred_rf = pred_RF_val,  y = y[val_samples, ]$beta,
                         sample = val_samples, type = "val")
  rbind(test_res, val_res)
}) -> pred_rf
pred_rf_res <- rbindlist(pred_rf)
cat(sprintf("RF test RMSE: %.4f\n",
            sqrt(mean((pred_rf_res[type == "test"]$pred_rf - pred_rf_res[type == "test"]$y)^2))))
saveRDS(pred_rf_res, sprintf("%s/pred_rf_res_chr%s.rds", res_dir, chr))

############################### xgboost ########################################
library(xgboost)
cat("xgboost\n")
pbmcapply::pbmclapply(EPIC$ID, function(x){
  generate_data(x, EPIC, HM450, predictors) -> tmp
  X <- tmp[[1]]; y <- tmp[[2]]
  xgb <- xgboost(data = X[train_samples, ], label = y[train_samples, ]$beta,
                 nrounds = 10, nthread = 1, objective = "reg:squarederror", verbose = 0)
  pred_XGB_test <- predict(xgb, X[test_samples, ])
  pred_XGB_val  <- predict(xgb, X[val_samples, ])
  test_res <- data.frame(probe = x, pred_xgb = pred_XGB_test, y = y[test_samples, ]$beta,
                         sample = test_samples, type = "test")
  val_res  <- data.frame(probe = x, pred_xgb = pred_XGB_val,  y = y[val_samples, ]$beta,
                         sample = val_samples, type = "val")
  rbind(test_res, val_res)
}) -> pred_xgb_test
pred_xgb_res <- rbindlist(pred_xgb_test)
cat(sprintf("XGB test RMSE: %.4f\n",
            sqrt(mean((pred_xgb_res[type == "test"]$pred_xgb - pred_xgb_res[type == "test"]$y)^2))))
saveRDS(pred_xgb_res, sprintf("%s/pred_xgb_res_chr%s.rds", res_dir, chr))

############################# KNN ##############################################
library(class)
cat("KNN\n")
pbmcapply::pbmclapply(EPIC$ID, function(x){
  generate_data(x, EPIC, HM450, predictors) -> tmp
  X <- tmp[[1]]; y <- tmp[[2]]
  knn(train = X[train_samples, ], test = X[c(val_samples, test_samples), ],
      cl = y[train_samples, ]$beta, k = 50, prob = T) -> knn
  pred_KNN <- as.numeric(as.matrix(knn))
  data.frame(probe = x,
             pred_knn = pred_KNN,
             y = c(y[val_samples, ]$beta, y[test_samples, ]$beta),
             type = c(rep("val", length(val_samples)), rep("test", length(test_samples))),
             sample = c(val_samples, test_samples))
}) -> pred_knn
pred_knn_res <- rbindlist(pred_knn)
cat(sprintf("KNN test RMSE: %.4f\n",
            sqrt(mean((pred_knn_res[type == "test"]$pred_knn - pred_knn_res[type == "test"]$y)^2))))
saveRDS(pred_knn_res, sprintf("%s/pred_knn_res_chr%s.rds", res_dir, chr))

############################# logistic regression #############################
cat("logistic regression\n")
pbmcapply::pbmclapply(EPIC$ID, function(x){
  generate_data(x, EPIC, HM450, predictors) -> tmp
  X <- tmp[[1]] %>% as.data.frame(); y <- tmp[[2]]
  meth <- as.vector(ifelse(y$beta > 0.5, 1, 0))
  X$meth <- meth
  glm(meth ~ ., data = X[train_samples, ], family = "binomial") -> glm
  pred_GLM_test <- predict(glm, X[test_samples, ], type = "response")
  pred_GLM_val  <- predict(glm, X[val_samples, ], type = "response")
  test_res <- data.frame(probe = x, pred_glm = pred_GLM_test, y = y[test_samples, ]$beta,
                         sample = test_samples, type = "test")
  val_res  <- data.frame(probe = x, pred_glm = pred_GLM_val,  y = y[val_samples, ]$beta,
                         sample = val_samples, type = "val")
  rbind(test_res, val_res)
}) -> pred_glm
pred_glm_res <- rbindlist(pred_glm)
cat(sprintf("GLM test RMSE: %.4f\n",
            sqrt(mean((pred_glm_res[type == "test"]$pred_glm - pred_glm_res[type == "test"]$y)^2))))
saveRDS(pred_glm_res, sprintf("%s/pred_glm_res_chr%s.rds", res_dir, chr))

############################# PFR ##############################################
if (Sys.getenv("RUN_PFR", "0") == "1") {
  cat("PFR (transfer)\n")
  source("src/CUE/Transfer_PFR.R")
} else {
  cat("PFR skipped (set RUN_PFR=1 to enable)\n")
}

cat("done\n")
