# implementation of CUE,
# ref: https://github.com/GangLiTarheel/CUE
library(data.table)
library(tidyverse)
options(mc.cores = 4)
options(ignore.interactive = TRUE)

# load the training data, chr
args = commandArgs(trailingOnly=TRUE)
# if (length(args) == 1) {
#   chr <- as.numeric(args[1])
# } else{
#   chr <- 20
# }
# print(chr)
# quit(save = "no")
chr <- 8
data_dir <- "tmp/AliciaKSmith_450K_EPIC/"
res_dir <- "res/AliciaKSmith_450K_EPIC_baseline/"

chr <- as.numeric(args[1])
data_dir <- args[2]
res_dir <- args[3]

HM450 <- fread(sprintf("%s/HM450_chr%d.csv", data_dir, chr))
EPIC <- fread(sprintf("%s/EPIC_chr%d.csv", data_dir, chr))

HM450 %>%
  drop_na() %>%
  arrange(MAPINFO) -> HM450

EPIC %>%
  drop_na() %>%
  arrange(MAPINFO) -> EPIC

# for EPIC data, use only the EPIC-only probes
EPIC %>%
  filter(ID %in% setdiff(EPIC$ID, HM450$ID)) -> EPIC
# sum(EPIC$ID %in% HM450$ID)

train_samples <- read.csv(sprintf("%s/train_samples.txt", data_dir), header = F)$V1
val_samples <- read.csv(sprintf("%s/val_samples.txt", data_dir), header = F)$V1
test_samples <- read.csv(sprintf("%s/test_samples.txt", data_dir), header = F)$V1

# for each probe in EPIC, find the upstream and downstream 25 probes in HM450,
# use the 25 probes to impute the EPIC probe

# get ID and MAPINFO for EPIC probes
EPIC %>% arrange(MAPINFO) -> EPIC
HM450 %>% arrange(MAPINFO) -> HM450

y_probes <- EPIC[, .(ID, MAPINFO)]
pbmcapply::pbmclapply(1:nrow(y_probes), function(i){
  # print(i)
  y_probes[i, .(ID, MAPINFO)] -> tmp
  # get the upstream and downstream 25 probes in HM450
  # based on MAPINFO
  # find the nearest probe in HM450 first
  dist <- abs(tmp$MAPINFO - HM450$MAPINFO)
  which.min(dist) -> idx
  if (idx <= 24) {
    return(data.frame(matrix(rep("X", 50), nrow = 1)))
  }
  if (idx > (nrow(HM450) - 25)) {
    return(data.frame(matrix(rep("X", 50), nrow = 1)))
  }
  # get the upstream and downstream 25 probes
  HM450[(idx - 24):(idx + 25), .(ID, MAPINFO)] -> tmp2
  matrix(tmp2$ID, nrow = 1) %>% as.data.frame()
}) -> near_probes

rbindlist(near_probes) -> X_probes
predictors <- cbind(y_probes, X_probes)

# remove invalid probes
predictors[predictors$X1 != "X"] -> predictors
setkey(predictors, ID)

EPIC %>%
  filter(ID %in% predictors$ID) -> EPIC
setkey(EPIC, ID)
EPIC <- EPIC[predictors$ID]
all(EPIC$ID == predictors$ID)

# for each probe, generate the X and y
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
  
  HM450[tmp_predictors,] %>%
    select(-chr, -MAPINFO) %>%
    select(ID, all_of(y$sample_name)) %>%
    as.data.frame() -> X
  rownames(X) <- X$ID
  X[,-1] -> X
  X <- t(X) %>% as.matrix()
  return(list(X, y))
}

############################# Random Forest ####################################
library(randomForest)
cat("Random Forest\n")
pbmcapply::pbmclapply(EPIC$ID, function(x){
  generate_data(x, EPIC, HM450, predictors) -> tmp
  X <- tmp[[1]]
  y <- tmp[[2]]
  # train
  randomForest(x = X[train_samples,], y = y[train_samples,]$beta,
               # xtest = X[val_samples,], ytest = y[val_samples,]$beta,
               ntree = 2,
               mtry = 3, verbose = 0) -> rf
  pred_RF_test<-predict(rf,X[test_samples,])
  pred_RF_val<-predict(rf,X[val_samples,])
  data.frame(probe = x, pred_rf = pred_RF_test, y = y[test_samples,]$beta,
             sample = test_samples) -> test_res
  data.frame(probe = x, pred_rf = pred_RF_val, y = y[val_samples,]$beta,
             sample = val_samples) -> val_res
  test_res$type <- "test"
  val_res$type <- "val"
  rbind(test_res, val_res)
}) -> pred_rf
pred_rf_res <- rbindlist(pred_rf)
# RMSE
sqrt(mean((pred_rf_res$pred_rf - pred_rf_res$y)^2))
saveRDS(pred_rf_res,
        sprintf("%s/pred_rf_res_chr%s.rds", res_dir, chr))
# quit(save = "no")

############################### xgboost ########################################
library(xgboost)
cat("xgboost\n")
pbmcapply::pbmclapply(EPIC$ID, function(x){
  # x <- "cg00007932"
  # print(x)
  generate_data(x, EPIC, HM450, predictors) -> tmp
  X <- tmp[[1]]
  y <- tmp[[2]]
  # train
  xgb <- xgboost(data = X[train_samples,], label = y[train_samples,]$beta,
                 nrounds = 10, nthread = 1,
                 objective = "reg:squarederror",verbose = 0)
  pred_XGB_test <- predict(xgb,X[test_samples,])
  data.frame(probe = x, pred_xgb = pred_XGB_test, y = y[test_samples,]$beta,
             sample = test_samples) -> test_res
  pred_XGB_val <- predict(xgb,X[val_samples,])
  data.frame(probe = x, pred_xgb = pred_XGB_val, y = y[val_samples,]$beta,
             sample = val_samples) -> val_res
  test_res$type <- "test"
  val_res$type <- "val"
  rbind(test_res, val_res)
}) -> pred_xgb_test
pred_xgb_res <- rbindlist(pred_xgb_test)
# # RMSE
sqrt(mean((pred_xgb_res$pred_xgb - pred_xgb_res$y)^2))
saveRDS(pred_xgb_res,
        sprintf("%s/pred_xgb_res_chr%s.rds", res_dir, chr))


############################# KNN ######################################
library(class)
cat("KNN\n")
pbmcapply::pbmclapply(EPIC$ID, function(x){
  generate_data(x, EPIC, HM450, predictors) -> tmp
  X <- tmp[[1]]
  y <- tmp[[2]]
  # train
  knn(train = X[train_samples,], test = X[c(val_samples,test_samples),],
      cl = y[train_samples,]$beta, k = 50, prob = T) -> knn
  pred_KNN_test <- as.numeric(as.matrix(knn))
  
  data.frame(probe = x, 
             pred_knn = pred_KNN_test, 
             y = c(y[val_samples,]$beta, y[test_samples,]$beta),
             type = c(rep("val", length(val_samples)), rep("test", length(test_samples))),
             sample = c(val_samples,test_samples))
}) -> pred_knn
pred_knn_res <- rbindlist(pred_knn)
# RMSE
sqrt(mean((pred_knn_res$pred_knn - pred_knn_res$y)^2))
saveRDS(pred_knn_res,
        sprintf("%s/pred_knn_res_chr%s.rds", res_dir, chr))

############################# logistic regression ######################
cat("logistic regression\n")
pbmcapply::pbmclapply(EPIC$ID, function(x){
  generate_data(x, EPIC, HM450, predictors) -> tmp
  X <- tmp[[1]] %>% as.data.frame()
  y <- tmp[[2]]
  meth <- ifelse(y$beta > 0.5, 1, 0)
  meth <- as.vector(meth)
  # glm_data <- X
  X$meth <- meth
  glm(meth ~ ., data = X[train_samples,], family = "binomial") -> glm
  pred_GLM_test <- predict(glm, X[test_samples,], type = "response")
  pred_GLM_val <- predict(glm, X[val_samples,], type = "response")
  data.frame(probe = x, pred_glm = pred_GLM_test, y = y[test_samples,]$beta,
             sample = test_samples) -> test_res
  data.frame(probe = x, pred_glm = pred_GLM_val, y = y[val_samples,]$beta,
             sample = val_samples) -> val_res
  test_res$type <- "test"
  val_res$type <- "val"
  rbind(test_res, val_res)
}) -> pred_glm
pred_glm_res <- rbindlist(pred_glm)
# RMSE
sqrt(mean((pred_glm_res$pred_glm - pred_glm_res$y)^2))
saveRDS(pred_glm_res,
        sprintf("%s/pred_glm_res_chr%s.rds", res_dir, chr))

############################# PFR ######################################
cat("PFR\n")
source("src/CUE/PFR.R")

# K = 4
# pred_CUE_equal_weighs <- 1/ K * (pred_Logistic_test + pred_KNN_test + pred_RF_test + pred_PFR_test)



