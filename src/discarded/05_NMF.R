library(NMF)
library(tidyverse)

encoding_450k <- read.csv("tmp/encoded/450k_chr20_4000.csv", row.names = 1)
encoding_EPIC <- read.csv("tmp/encoded/EPIC_chr20_8000.csv", row.names = 1)

# Test data ['AFB182' 'EUB146' 'AFB155' 'EUB046' 'EUB094' 'AFB035' 'EUB134' 'AFB022'
# 'AFB082' 'EUB153' 'EUB120' 'AFB048' 'EUB043' 'EUB003' 'EUB052' 'AFB053'
# 'AFB163']
test_samples <- c("AFB182", "EUB146", "AFB155", "EUB046", "EUB094", "AFB035",
                  "EUB134", "AFB022", "AFB082", "EUB153", "EUB120", "AFB048",
                  "EUB043", "EUB003", "EUB052", "AFB053", "AFB163")

encoding_450k_test <- encoding_450k[test_samples, ]
encoding_EPIC_test <- encoding_EPIC[test_samples, ]

# remove test samples, rownames are sample names
encoding_450k <- encoding_450k[setdiff(rownames(encoding_450k), test_samples), ]
encoding_EPIC <- encoding_EPIC[setdiff(rownames(encoding_EPIC), test_samples), ]

# NMF
set.seed(123)
# make sure the number of samples are the same
rownames(encoding_450k) == rownames(encoding_EPIC)
mat <- cbind(encoding_450k, encoding_EPIC) %>% as.matrix()
dim(mat)

res <- nmf(t(mat), rank = 50)
H <- res@fit@H
W <- res@fit@W
# mutiply H and W to get the reconstructed matrix
reconstructed_mat <- W %*% H
dim(reconstructed_mat)
dim(mat)
error <- t(mat) - reconstructed_mat
sum(abs(error)) / nrow(mat) / ncol(mat)

# 450K -> EPIC, encoding_450k_test <- W1 %*% H0, to get H0
W1 <- W[1:ncol(encoding_450k), ]
W2 <- W[(ncol(encoding_450k) + 1):nrow(W), ]
dim(encoding_450k_test)
dim(W1)
dim(W2)

res2 <- nmf(t(encoding_450k_test), W = W1, rank = 50)
res2@fit@H %>% dim()
H0 <- res2@fit@H
# H0 is the encoding of 450K test samples in EPIC space
H0 %>% dim()

W2 %>% dim()

W2 %*% H0 %>% t() -> reconstructed_mat_450k_to_EPIC

abs(reconstructed_mat_450k_to_EPIC - encoding_EPIC_test) %>% sum() / nrow(encoding_EPIC_test) / ncol(encoding_EPIC_test)



