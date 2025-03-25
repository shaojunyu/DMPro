library(tidyverse)
library(data.table)
library(ggvenn)
library(ggpubr)
library(ggsci)

# train log
hm450_train <- fread("tmp/450k_train_log.txt", header = F)
hm450_train$V2 %>% str_remove("Chr20 Epoch \\[") %>%
  str_remove("/300]") %>%
  as.numeric() -> hm450_train$epoch
hm450_train$V3 %>% str_remove("Train Loss: ") %>%
  as.numeric() -> hm450_train$train_loss
hm450_train$V4 %>% str_remove("Val Loss: ") %>%
  as.numeric() -> hm450_train$val_loss

hm450_train %>%
  pivot_longer(cols = c("train_loss", "val_loss"), names_to = "type", values_to = "loss") %>%
  ggplot(aes(x = epoch, y = loss, color = type)) +
  geom_line(linewidth = 1) +
  xlab("Training Epoch") +
  ylab("Loss (MAE)") + ggtitle("HM450 Training") +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"))

epic_train <- fread("tmp/epic_train_log.txt", header = F)
epic_train$V2 %>% str_remove("Chr20 Epoch \\[") %>%
  str_remove("/300]") %>%
  as.numeric() -> epic_train$epoch
epic_train$V3 %>% str_remove("Train Loss: ") %>%
  as.numeric() -> epic_train$train_loss
epic_train$V4 %>% str_remove("Val Loss: ") %>%
  as.numeric() -> epic_train$val_loss

epic_train %>%
  pivot_longer(cols = c("train_loss", "val_loss"), names_to = "type", values_to = "loss") %>%
  ggplot(aes(x = epoch, y = loss, color = type)) +
  geom_line(linewidth = 1) +
  xlab("Training Epoch") +
  ylab("Loss (MAE)") + ggtitle("EPIC Training") +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"))


pfr_res <- readRDS("res/CUE/pred_pfr_res.rds")
glm_res <- readRDS("res/CUE/pred_glm_res.rds")
knn_res <- readRDS("res/CUE/pred_knn_res.rds")
rf_res <- readRDS("res/CUE/pred_rf_res.rds")
xgb_res <- readRDS("res/CUE/pred_xgb_res.rds")

train_samples <- read.csv("res/NMF/train_samples.txt", header = F)$V1
test_samples <- read.csv("tmp/NMF/test_samples.txt", header = F)$V1
val_samples <- read.csv("tmp/NMF/val_samples.txt", header = F)$V1

# using chr20
epic_20 <- fread("tmp/NMF/EPIC_chr20.csv") %>% drop_na()
hm450_20 <- fread("tmp/NMF/HM450_chr20.csv") %>% drop_na()
setkey(epic_20, ID)
setkey(hm450_20, ID)

# get the epic-specific probes
setdiff(epic_20$ID, hm450_20$ID) -> epic_specific
shared_probes <- intersect(epic_20$ID, hm450_20$ID)

decoded <- fread("tmp/decoded/EPIC_chr20_8000_bk.csv") %>% drop_na()
setkey(decoded, probe)
decoded <- decoded[epic_20$ID]

lapply(test_samples, function(x){
  # RMSE for epic-specific probes
  decoded[epic_specific, ..x] -> tmp1
  epic_20[epic_specific, ..x] -> tmp2
  rmse1 = sqrt(mean((tmp1[[x]] - tmp2[[x]])^2))
  MAE1 = mean(abs(tmp1[[x]] - tmp2[[x]]))
  # sprintf("RMSE for %s is %f", x, rmse)
  
  # RMSE for shared probes
  hm450_20[shared_probes, ..x] -> tmp1
  epic_20[shared_probes, ..x] -> tmp2
  rmse2 = sqrt(mean((tmp1[[x]] - tmp2[[x]])^2))
  MAE2 = mean(abs(tmp1[[x]] - tmp2[[x]]))
  
  # RMSE for shared probes
  decoded[shared_probes, ..x] -> tmp1
  epic_20[shared_probes, ..x] -> tmp2
  rmse3 = sqrt(mean((tmp1[[x]] - tmp2[[x]])^2))
  
  # overall RMSE
  decoded[epic_20$ID, ..x] -> tmp1
  epic_20[epic_20$ID, ..x] -> tmp2
  rmse4 = sqrt(mean((tmp1[[x]] - tmp2[[x]])^2))
  # sprintf("%s, Epic-RMSE %f, shared RMSE %f, %f", x, rmse1, rmse2, rmse3)
  data.frame(sample = x, Our_rmse = rmse1, Shared_rmse = rmse2)
  # data.frame(sample = x, Our_MAE = MAE1, Shared_MAE = MAE2)
}) %>% rbindlist() -> our_rmse

# RMSE
knn_res %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(KNN_rmse = sqrt(mean((y - pred_knn)^2))) -> knn_rmse

rf_res %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(rmse = sqrt(mean((y - pred_rf)^2))) %>%
  summarise(mean(rmse))

pfr_res %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(PFR_rmse = sqrt(mean((y - pred_PFR)^2))) %>%
  drop_na() -> pfr_rmse
  summarise(mean(rmse))

xgb_res %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(XGBoost_rmse = sqrt(mean((y - pred_xgb)^2))) -> xgb_rmse
  summarise(mean(rmse))

glm_res %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(rmse = sqrt(mean((y - pred_glm)^2))) %>%
  summarise(mean(rmse))

# CUE
1 / 2 * (knn_res$pred_knn + xgb_res$pred_xgb) -> pred_cue
xgb_res$pred_cue <- pred_cue
xgb_res %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(CUE_rmse = sqrt(mean((y - pred_cue)^2))) -> cue_rmse

# res_df <- data.frame(
#   method = c("Our", "KNN", "XGBoost", "PFR", "CUE"),
#   rmse = c(0.03918268, 0.0518, 0.0435, 0.0512, 0.0430)
# )
# combine all the rmse results
res_df <- knn_rmse %>%
  left_join(pfr_rmse) %>%
  left_join(xgb_rmse) %>%
  left_join(cue_rmse) %>%
  left_join(our_rmse) %>%
  pivot_longer(cols = -sample, names_to = "method", values_to = "rmse") %>%
  mutate(method = str_replace(method, "_rmse", ""))

# violin plot
unique(res_df$method)
res_df %>%
  mutate(method = factor(method, levels = c(
    "KNN", "PFR", "XGBoost", "CUE", "Our", "Shared"
  ))) %>%
  ggplot(aes(x = method, y = rmse)) +
  geom_boxplot(color = "black",
               fill = "white", width = 0.5) +
  theme_classic() + xlab("") + ylab("RMSE") +
  # scale_y_continuous(limits = c(0, 0.055)) +
  # scale_color_nejm() +
  theme(axis.text.x = element_text(color = "black", face = "bold",
                                   size = 12),
        axis.text.y = element_text(size = 12, color = "black")) -> p1

res_df %>%
  drop_na() %>%
  group_by(method) %>%
  summarise(rmse = mean(rmse)) %>%
  arrange(desc(rmse)) %>%
  mutate(method = factor(method, levels = method)) %>%
  ggplot(aes(x = method, y = rmse)) +
  geom_bar(fill = "steelblue",
           stat = "identity", 
           width = 0.6) +
  geom_text(aes(label = round(rmse, 4)), color = "black", 
            fontface = "bold",
            size = 4.5,
            vjust = -0.2) +
  geom_hline(yintercept = mean(our_rmse$Our_rmse), linetype = "dashed", color = "red") +
  theme_bw() + xlab("") + ylab("RMSE") +
  # scale_fill_nejm() +
  # scale_y_continuous(limits = c(0, 0.056)) +
  theme(axis.text.x = element_text(color = "black", face = "bold",
                                   size = 12),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none") -> p2
library(patchwork)
p1 | p2 + plot_annotation(title = )
p2 + ggtitle("Chr_19") +
  theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold", size = 14))



# library(ggpubr)
# ggbarplot(data = res_df,
#                   x = "method", y = "rmse", fill = "method")

# check the rownames
# decoded$probe == original$ID


data.frame(original = original$AFB045,
           decoded = decoded$AFB045,
           probe = decoded$probe) %>%
  na.omit() %>%
  ggplot(aes(x = original, y = decoded)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

# calculate the MES by using the epic-specific probes
data.frame(original = original$AFB045,
           decoded = decoded$AFB045,
           probe = decoded$probe) %>%
  filter(probe %in% epic_specific)

epic_20[,..test_samples] %>%
  as.matrix() -> original_mat

decoded[,..test_samples] %>%
  as.matrix() -> decoded_mat

rownames(decoded_mat) <- decoded$probe
rownames(original_mat) <- epic_20$ID
rownames(original_mat) == rownames(decoded_mat)
colnames(decoded_mat) <- paste0("decoded_", colnames(decoded_mat))
colnames(original_mat) <- paste0("original_", colnames(original_mat))

mat <- cbind(original_mat, decoded_mat)
as.data.frame(mat) %>% drop_na() -> mat
mat[epic_specific,] -> mat
as.data.frame(mat) %>% drop_na() -> mat

mat[epic_specific,] -> mat
as.data.frame(mat) %>% drop_na() -> mat

# z-score normalization for each column
mat %>% scale() -> mat_scaled

# mat$fake <- mean(mat$original_AFB022)
# na.omit(mat)
# mat <- matrix(mat)
pairwise_MAE <- function(mat) {
  n <- ncol(mat)
  res <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      res[i, j] <- mean(abs(mat[,i] - mat[,j]))
    }
  }
  res
  colnames(res) <- colnames(mat)
  rownames(res) <- colnames(mat)
  res
}
pairwise_MAE(mat = mat) -> MAE
pheatmap::pheatmap(MAE, display_numbers = T,
                   cluster_rows = F, cluster_cols = F)

# MSE
pairwise_MSE <- function(mat) {
  n <- ncol(mat)
  res <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      res[i, j] <- mean((mat[,i] - mat[,j])^2)
    }
  }
  res
  colnames(res) <- colnames(mat)
  rownames(res) <- colnames(mat)
  res
}
pairwise_MSE(mat = mat) -> MSE
pheatmap::pheatmap(MSE, display_numbers = T,
                   cluster_rows = F, cluster_cols = F)

# RMSE
pairwise_RMSE <- function(mat) {
  n <- ncol(mat)
  res <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      res[i, j] <- sqrt(mean((mat[,i] - mat[,j])^2))
    }
  }
  res
  colnames(res) <- colnames(mat)
  rownames(res) <- colnames(mat)
  res
}
pairwise_RMSE(mat = mat) -> RMSE
pheatmap::pheatmap(RMSE, display_numbers = T,
                   cluster_rows = F, cluster_cols = F)

cor_mat <- cor(mat_scaled)
pheatmap::pheatmap(cor_mat, display_numbers = T,
                   cluster_rows = F, cluster_cols = F)
# pairwise_MSE(mat = mat) -> original_MSE
# pheatmap::pheatmap(original_MSE, display_numbers = T,
#                    cluster_rows = F, cluster_cols = F)




# data.frame(original = original$AFB022,
#            decoded = decoded$AFB022,
#            probe = decoded$probe) %>%
#   filter(probe %in% epic_specific) %>%
#   na.omit() %>%
#   ggplot(aes(x = original, y = decoded)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   stat_cor(method = "pearson") +
#   theme_bw() +
#   labs(title = "EPIC-specific probes") +
#   xlab("Original EPIC") +
#   ylab("Decoded from 450K")



