n = 200              # number of sampels
x <- matrix(rbeta(n*50,1,1),ncol=n) # upperstream 25 probes and downstream 25 probes measured by HM450
y <- rbeta(n,1,1)    # the target HM850 probe

Meth<-data.frame(cbind(y,t(x)))
train=1:160
test=161:200

source("tmp/CUE/R/refund_lib.R")
load("tmp/CUE/PTSD/Annotations.RData")
#load("Data/Annotations.RData")
p=dim(annotation.450K)[1]
X_450 = matrix(rbeta(n*p,1,1),ncol=n) 
rownames(X_450) <- rownames(annotation.450K)
X_450 = log2(X_450/(1-X_450))

Y=log2(y/(1-y))
X=log2(x/(1-x))

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
  test.dens[[j]] <- X_450[annotation.450K[, "Relation_to_UCSC_CpG_Island"] == l,test]
  test.funcs[[j]] <- apply(test.dens[[j]], 2, function(x) {density(x,from=-13.38757,to=13.38757)$y})
  train.dens[[j]] <- X_450[annotation.450K[, "Relation_to_UCSC_CpG_Island"] == l,train] 
  train.funcs[[j]] <- apply(train.dens[[j]], 2, function(x) {density(x,from=-13.38757,to=13.38757)$y})   
}

## random select a target probe Y to impute:
# target_probe = rownames(annotation.EPIC)[sample.int(dim(annotation.EPIC),1)]
target_probe = rownames(annotation.EPIC)[1]
j <- paste(annotation.EPIC[target_probe,"Relation_to_UCSC_CpG_Island"])
if (j == "") {
  j <- "NA"
}
fit = new_pfr(unlist(Y[train]),
              t(X[,train]), 
              t(train.funcs[[j]]))

X.test <- create_predictors(t(X[,test]), t(test.funcs[[j]]))
Y.test <- X.test %*% fit$coefs
pred_PFR_test <- 2 ^ Y.test / (2 ^ Y.test + 1)
