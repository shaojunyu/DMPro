# implementatio of the PFR (penalized functional regression) method
# ref: https://github.com/Leonardo0628/pfr
# git clone https://github.com/Leonardo0628/pfr to tmp/pfr
# git clone https://github.com/GangLiTarheel/CUE to tmp/CUE
# This script is adapted from the pfr.R script in the pfr repository.
# The original script is used to impute 450K data using 27K data.
# Here we use it to impute EPIC data using 450K data.
library(data.table)
library(tidyverse)
# source("tmp/pfr/refund_lib.R")

# Meth<-data.frame(cbind(y$beta,t(X)))
source("tmp/CUE/R/refund_lib.R")

if(imputation_type == "HM450_EPIC"){
  X_450_train <- fread(sprintf("%s/HM450_train.csv", data_dir))
  X_450_test <- fread(sprintf("%s/HM450_test.csv", data_dir))
  X_450_val <- fread(sprintf("%s/HM450_val.csv", data_dir))
}

if(imputation_type == "EPIC_V1_V2"){
  X_450_train <- fread(sprintf("%s/EPIC_V1_train.csv", data_dir))
  X_450_test <- fread(sprintf("%s/EPIC_V1_test.csv", data_dir))
  X_450_val <- fread(sprintf("%s/EPIC_V1_val.csv", data_dir))
}


ID = X_450_train$ID
all(X_450_train$ID == X_450_test$ID)
X_450_train <- X_450_train %>%
  select(-chr, -MAPINFO, -ID) -> X_450_train
X_450_test <- X_450_test %>%
  select(-chr, -MAPINFO, -ID) -> X_450_test
X_450_val <- X_450_val %>%
  select(-chr, -MAPINFO, -ID) -> X_450_val
X_450 <- cbind(X_450_train, X_450_test, X_450_val) %>% as.matrix()
train_idx <- which(colnames(X_450) %in% train_samples)
test_samples <- c(test_samples, val_samples)
test_idx <- which(colnames(X_450) %in% test_samples)

if(imputation_type == "HM450_EPIC"){
  annotation.450K <- fread(cmd = 'grep -v ^# data/annotations/GPL13534\\(450K\\)-11288.txt') %>%
    as.data.frame()
  rownames(annotation.450K) <- annotation.450K$ID
  annotation.450K <- annotation.450K[ID,]
}

if(imputation_type == "EPIC_V1_V2"){
  annotation.450K <- fread(cmd = 'grep -v ^# data/annotations/GPL33022\\(EPICv2\\)-233683.txt') %>%
    as.data.frame()
  rownames(annotation.450K) <- annotation.450K$ID
  annotation.450K <- annotation.450K[ID,]
}

rownames(X_450) <- rownames(annotation.450K)
X_450 = log2(X_450/(1-X_450))

isl_group <- c("", "Island", "N_Shelf", "N_Shore", "S_Shore", "S_Shelf")

train.dens <- list()
train.funcs <- list()
test.dens <- list()
test.funcs <- list()

for (l in isl_group) {
  if (l == "") {
    j <- "NA"
  } else {
    j <- l
  }
  # max(X)==13.28757
  test.dens[[j]] <- X_450[annotation.450K[, "Relation_to_UCSC_CpG_Island"] == l,test_idx]
  test.funcs[[j]] <- apply(test.dens[[j]], 2, function(x) {density(x,from=-13.38757,to=13.38757)$y})
  train.dens[[j]] <- X_450[annotation.450K[, "Relation_to_UCSC_CpG_Island"] == l,train_idx] 
  train.funcs[[j]] <- apply(train.dens[[j]], 2, function(x) {density(x,from=-13.38757,to=13.38757)$y})   
}

# data/annotations/GPL23976(850K)-56952.txt
annotation.EPIC <- fread(cmd = 'grep -v ^# data/annotations/GPL33022\\(EPICv2\\)-233683.txt') %>%
  as.data.frame()
annotation.EPIC <- annotation.EPIC[!duplicated(annotation.EPIC$Name),]
row.names(annotation.EPIC) <- annotation.EPIC$Name
annotation.EPIC <- annotation.EPIC[predictors$ID,]

lapply(EPIC$ID, function(target_probe){
  # print(target_probe)
  generate_data(target_probe, EPIC, HM450, predictors) -> tmp
  x <- tmp[[1]]
  y <- tmp[[2]]
  Y=log2(y$beta/(1-y$beta))
  # prevent 0 and 1
  x[x == 0] <- 0.000001
  x[x == 1] <- 0.9999999
  X=t(log2(x/(1-x)))
  j <- paste(annotation.EPIC[target_probe, "Relation_to_UCSC_CpG_Island"])
  if (j == "") {
    j <- "NA"
  }
  j <- str_split(j,";")[[1]][1]
  
  fit <- try( new_pfr(unlist(Y[train_idx]),
                t(X[,train_idx]),
                t(train.funcs[[j]])), TRUE)
  if (class(fit) == "try-error") {
    NA
  }else{
    X.test <- create_predictors(t(X[,test_idx]), t(test.funcs[[j]]))
    Y.test <- X.test %*% fit$coefs
    pred_PFR_test <- 2 ^ Y.test / (2 ^ Y.test + 1)
    data.frame(probe = target_probe, pred_PFR = pred_PFR_test, y = y[test_samples,]$beta,
                 sample = test_samples)
  }
}) -> pred_pfr
sum(!is.na(pred_pfr))

# remve NA
pred_pfr <- pred_pfr[!is.na(pred_pfr)]
pred_pfr_res <- rbindlist(pred_pfr[!is.na(pred_pfr)])

#RMSE
sqrt(mean((pred_pfr_res$pred_PFR - pred_pfr_res$y)^2))
saveRDS(pred_pfr_res,
        sprintf("%s/pred_pfr_res_chr%s.rds", res_dir, chr))
