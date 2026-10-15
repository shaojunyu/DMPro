# PFR for the cross-dataset (transfer) setting: fit on the SOURCE dataset's
# train samples, predict the TARGET dataset's test samples (+ source val, for CUE).
# Sourced by src/CUE/Transfer_baseline.R when RUN_PFR=1.
# Samples are indexed by name, not position, since they come from two datasets.
# ref: https://github.com/Leonardo0628/pfr , https://github.com/GangLiTarheel/CUE

library(mgcv)
source("tmp/CUE/R/refund_lib.R")

pfr_pred_samples <- c(test_samples, val_samples)   # target test + source val

# ---- genome-wide HM450 matrix spanning all samples we need -----------------
read_split <- function(path, want) {
  hdr <- names(fread(path, nrows = 0))
  keep <- intersect(want, hdr)
  if (length(keep) == 0) return(NULL)
  fread(path, select = c("ID", keep))
}

x450_src_train <- read_split(sprintf("%s/HM450_train.csv", source_dir), train_samples)
x450_src_val   <- read_split(sprintf("%s/HM450_val.csv",   source_dir), val_samples)
x450_tgt_test  <- read_split(sprintf("%s/HM450_test.csv",  target_dir), test_samples)
stopifnot(!is.null(x450_src_train), !is.null(x450_src_val), !is.null(x450_tgt_test))

# common genome-wide probes across the three sources, in a single fixed order
common_450 <- Reduce(intersect, list(x450_src_train$ID, x450_src_val$ID, x450_tgt_test$ID))
stopifnot(length(common_450) > 0)
align <- function(d) {
  setkey(d, ID)
  d <- d[common_450]
  as.matrix(d[, -1, with = FALSE])
}
X_450 <- cbind(align(x450_src_train), align(x450_src_val), align(x450_tgt_test))
rownames(X_450) <- common_450
rm(x450_src_train, x450_src_val, x450_tgt_test)

train_cols <- intersect(train_samples, colnames(X_450))
pred_cols  <- intersect(pfr_pred_samples, colnames(X_450))
cat(sprintf("PFR chr%s: genome-wide 450K probes=%d  train(src)=%d  pred(tgt test+src val)=%d\n",
            chr, nrow(X_450), length(train_cols), length(pred_cols)))
stopifnot(length(train_cols) > 0, length(pred_cols) > 0)

# logit transform; clamp so exact 0/1 betas do not become +-Inf
X_450[X_450 <= 0] <- 1e-6
X_450[X_450 >= 1] <- 1 - 1e-7
X_450 <- log2(X_450 / (1 - X_450))

# ---- island annotation for the predictor (450K) probes ---------------------
annotation.450K <- fread(cmd = 'grep -v ^# data/annotations/GPL13534\\(450K\\)-11288.txt') %>%
  as.data.frame()
rownames(annotation.450K) <- annotation.450K$ID
annotation.450K <- annotation.450K[rownames(X_450), ]
isl_450 <- annotation.450K[, "Relation_to_UCSC_CpG_Island"]
isl_450[is.na(isl_450)] <- ""      # probes absent from the annotation -> open sea

isl_group <- c("", "Island", "N_Shelf", "N_Shore", "S_Shore", "S_Shelf")

train.funcs <- list()
test.funcs <- list()
for (l in isl_group) {
  j <- if (l == "") "NA" else l
  rows <- which(isl_450 == l)
  if (length(rows) < 2) next          # too few probes to estimate a density
  train.funcs[[j]] <- apply(X_450[rows, train_cols, drop = FALSE], 2,
                            function(x) density(x, from = -13.38757, to = 13.38757)$y)
  test.funcs[[j]]  <- apply(X_450[rows, pred_cols,  drop = FALSE], 2,
                            function(x) density(x, from = -13.38757, to = 13.38757)$y)
}

# ---- island annotation for the target (EPIC) probes -----------------------
# EPICv2 annotation: the EPIC (850K) one has no island column
annotation.EPIC <- fread(cmd = 'grep -v ^# data/annotations/GPL33022\\(EPICv2\\)-233683.txt') %>%
  as.data.frame()
annotation.EPIC <- annotation.EPIC[!duplicated(annotation.EPIC$Name), ]
rownames(annotation.EPIC) <- annotation.EPIC$Name
annotation.EPIC <- annotation.EPIC[predictors$ID, ]

# ---- per-probe functional regression --------------------------------------
lapply(EPIC$ID, function(target_probe) {
  tmp <- generate_data(target_probe, EPIC, HM450, predictors)
  x <- tmp[[1]]; y <- tmp[[2]]
  x[x <= 0] <- 1e-6
  x[x >= 1] <- 1 - 1e-7
  X <- t(log2(x / (1 - x)))                     # predictors x samples
  b <- y$beta
  b[b <= 0] <- 1e-6
  b[b >= 1] <- 1 - 1e-7
  Y <- setNames(log2(b / (1 - b)), y$sample_name)

  j <- paste(annotation.EPIC[target_probe, "Relation_to_UCSC_CpG_Island"])
  if (j == "" || j == "NA" || is.na(j)) j <- "NA"
  j <- str_split(j, ";")[[1]][1]
  if (!j %in% names(train.funcs)) j <- "NA"
  if (!j %in% names(train.funcs)) return(NA)

  tr <- intersect(train_cols, colnames(X))
  te <- intersect(pred_cols,  colnames(X))
  if (length(tr) == 0 || length(te) == 0) return(NA)

  fit <- try(new_pfr(unlist(Y[tr]),
                     t(X[, tr, drop = FALSE]),
                     t(train.funcs[[j]][, tr, drop = FALSE])), TRUE)
  if (inherits(fit, "try-error")) return(NA)

  X.test <- create_predictors(t(X[, te, drop = FALSE]),
                              t(test.funcs[[j]][, te, drop = FALSE]))
  Y.test <- X.test %*% fit$coefs
  pred_PFR_test <- 2^Y.test / (2^Y.test + 1)
  data.frame(probe = target_probe,
             pred_PFR = as.numeric(pred_PFR_test),
             y = y[te, ]$beta,
             sample = te)
}) -> pred_pfr

pred_pfr <- pred_pfr[!is.na(pred_pfr)]
cat(sprintf("PFR chr%s: %d/%d probes fitted\n", chr, length(pred_pfr), nrow(EPIC)))
pred_pfr_res <- rbindlist(pred_pfr)

cat(sprintf("PFR chr%s test+val RMSE: %.4f\n", chr,
            sqrt(mean((pred_pfr_res$pred_PFR - pred_pfr_res$y)^2))))
saveRDS(pred_pfr_res, sprintf("%s/pred_pfr_res_chr%s.rds", res_dir, chr))
