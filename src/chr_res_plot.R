library(tidyverse)
library(data.table)

# load the training data, chr
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1) {
  chr <- as.numeric(args[1])
} else{
  chr <- 20
}

pfr_res <- readRDS(sprintf("res/CUE/pred_pfr_res_chr%s.rds", chr))
glm_res <- readRDS(sprintf("res/CUE/pred_glm_res_chr%s.rds", chr))
knn_res <- readRDS(sprintf("res/CUE/pred_knn_res_chr%s.rds", chr))
rf_res <- readRDS(sprintf("res/CUE/pred_rf_res_chr%s.rds", chr))
xgb_res <- readRDS(sprintf("res/CUE/pred_xgb_res_chr%s.rds", chr))

train_samples <- read.csv("tmp/processed/train_samples.txt", header = F)$V1
test_samples <- read.csv("tmp/processed/test_samples.txt", header = F)$V1
val_samples <- read.csv("tmp/processed/val_samples.txt", header = F)$V1

# read data
EPIC <- fread(sprintf("tmp/processed/EPIC_chr%s.csv", chr)) %>% drop_na()
HM450 <- fread(sprintf("tmp/processed/HM450_chr%s.csv", chr)) %>% drop_na()
setkey(EPIC, ID)
setkey(HM450, ID)

# get the epic-specific probes
setdiff(EPIC$ID, HM450$ID) -> epic_specific
shared_probes <- intersect(EPIC$ID, HM450$ID)

decoded <- fread(sprintf("tmp/decoded/EPIC_chr%s_8000.csv", chr)) %>% drop_na()
setkey(decoded, probe)
decoded <- decoded[EPIC$ID]


lapply(test_samples, function(x){
  # x <- "AFB045"
  # RMSE for epic-specific probes
  decoded[epic_specific, ..x] -> tmp1
  EPIC[epic_specific, ..x] -> tmp2
  rmse1 = sqrt(mean((tmp1[[x]] - tmp2[[x]])^2))
  MAE1 = mean(abs(tmp1[[x]] - tmp2[[x]]))
  # sprintf("RMSE for %s is %f", x, rmse)
  
  # RMSE for shared probes
  HM450[shared_probes, ..x] -> tmp1
  EPIC[shared_probes, ..x] -> tmp2
  rmse2 = sqrt(mean((tmp1[[x]] - tmp2[[x]])^2))
  MAE2 = mean(abs(tmp1[[x]] - tmp2[[x]]))
  
  # RMSE for shared probes
  decoded[shared_probes, ..x] -> tmp1
  EPIC[shared_probes, ..x] -> tmp2
  rmse3 = sqrt(mean((tmp1[[x]] - tmp2[[x]])^2))
  
  # overall RMSE
  decoded[EPIC$ID, ..x] -> tmp1
  EPIC[EPIC$ID, ..x] -> tmp2
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
  summarise(RF_rmse = sqrt(mean((y - pred_rf)^2))) -> rf_rmse
  # summarise(mean(rmse)) -> rf_rmse

pfr_res %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(PFR_rmse = sqrt(mean((y - pred_PFR)^2))) %>%
  drop_na() -> pfr_rmse
# summarise(mean(rmse))

xgb_res %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(XGBoost_rmse = sqrt(mean((y - pred_xgb)^2))) -> xgb_rmse
# summarise(mean(rmse))

glm_res %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(rmse = sqrt(mean((y - pred_glm)^2))) %>%
  summarise(mean(rmse))


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
  select(probe, sample, pred_cue, y.x) -> cue_res

# cbind(cue_res$pred_knn, cue_res$pred_xgb, cue_res$pred_rf) %>%
#   rowMeans(na.rm = T) -> cue_res$pred_cue

cue_res %>%
  replace_na(list(pred_cue = 0)) %>%
  filter(probe %in% epic_specific,
         sample %in% test_samples) %>%
  group_by(sample) %>%
  summarise(CUE_rmse = sqrt(mean((y.x - pred_cue)^2))) -> cue_rmse



# combine all the rmse results
res_df <- knn_rmse %>%
  left_join(pfr_rmse) %>%
  left_join(xgb_rmse) %>%
  left_join(cue_rmse) %>%
  left_join(our_rmse) %>%
  left_join(rf_rmse) %>%
  pivot_longer(cols = -sample, names_to = "method", values_to = "rmse") %>%
  mutate(method = str_replace(method, "_rmse", ""))

# violin plot
unique(res_df$method)
res_df %>%
  drop_na() %>%
  group_by(method) %>%
  summarise(rmse = mean(rmse)) %>%
  arrange(desc(rmse)) %>%
  mutate(method = factor(method, levels = method)) -> rank
res_df %>%
  mutate(method = factor(method, levels = rank$method)) %>%
  ggplot(aes(x = method, y = rmse)) +
  geom_boxplot(color = "black",
               fill = "white", width = 0.5) +
  theme_classic() + xlab("") + ylab("RMSE") +
  # scale_y_continuous(limits = c(0, 0.055)) +
  # scale_color_nejm() +
  theme(axis.text.x = element_text(color = "black", face = "bold",
                                   size = 12),
        axis.text.y = element_text(size = 12, color = "black")) -> p1
p1
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
p1 | (p2 + theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold", size = 14)))
p1 + ggtitle(sprintf("Chr_%s", chr)) +
  theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold", size = 14))
export::graph2pdf(last_plot(),
                  sprintf("res/plot/chr%s_rmse_boxplot.pdf", chr),
                  width = 6, height = 5)
p2 + ggtitle(sprintf("Chr_%s", chr)) +
  theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold", size = 14))
export::graph2pdf(last_plot(),
                  sprintf("res/plot/chr%s_rmse_barplot.pdf", chr),
                  width = 6, height = 5)


# example scatter plot
sample_name = "EUB197"
df <- data.frame(
  observed = EPIC[epic_specific, ..sample_name] %>% unlist() %>% as.numeric(),
  imputed = decoded[epic_specific, ..sample_name] %>% unlist() %>% as.numeric()
)
# ggplot(df, aes(x = observed, y = imputed)) +
#   geom_point(color = "black", alpha = 0.5) +
#   # add diagonal line
#   geom_abline(intercept = 0, slope = 1, color = "red") +
#   # add pearson correlation 
#   geom_label(aes(x = 0.1, y = 0.9, label = sprintf("r = %.3f", cor(observed, imputed, method = "pearson"))),
#             color = "black", size = 4.5, fontface = "bold") +
#   theme_bw() +
#   theme(axis.title = element_text(color = "black", face = "bold",
#                                   size = 13)) +
#   ggtitle(sprintf("Sample:%s, Chr_%s",sample_name, chr)) +
  # theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold", size = 14))

library(ggpubr)
ggscatter(df, x = "observed", y = "imputed", color = "grey", alpha = 0.5,
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Observed", ylab = "Imputed") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", face = "bold",
                                   size = 12)) +
  ggtitle(sprintf("Sample:%s, Chr_%s",sample_name, chr)) +
  theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold", size = 14))
export::graph2png(last_plot(),
                  sprintf("res/plot/chr%s_scatter_%s.png", chr, sample_name),
                  width = 6, height = 5)

