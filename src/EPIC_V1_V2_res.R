library(tidyverse)
library(data.table)
library(ggpubr)
library(ggsci)

data_dir <- "tmp/EPIC_V1_V2/"
baseline_dir <- "res/EPIC_V1_V2_baseline/"
decoded_dir <- "res/EPIC_V1_V2_decoded/"

train_samples <- read.csv(paste0(data_dir, "/train_samples.txt"), header = F)$V1
test_samples <- read.csv(paste0(data_dir, "/test_samples.txt"), header = F)$V1
val_samples <- read.csv(paste0(data_dir, "/val_samples.txt"), header = F)$V1

# read data for each chr
options(mc.cores = 22)
pbmcapply::pbmclapply(20:20, function(chr){
  # chr <- 20
  print(chr)
  pfr_res <- readRDS(sprintf("%s/pred_pfr_res_chr%s.rds", baseline_dir, chr))
  glm_res <- readRDS(sprintf("%s/pred_glm_res_chr%s.rds", baseline_dir, chr))
  knn_res <- readRDS(sprintf("%s/pred_knn_res_chr%s.rds", baseline_dir, chr))
  rf_res <- readRDS(sprintf("%s/pred_rf_res_chr%s.rds", baseline_dir, chr))
  xgb_res <- readRDS(sprintf("%s/pred_xgb_res_chr%s.rds", baseline_dir, chr))
  cue_res <- readRDS(sprintf("%s/pred_cue_res_chr%s.rds",baseline_dir, chr))
  
  # read data
  EPIC <- fread(sprintf("%s/EPIC_V2_chr%s.csv", data_dir, chr)) %>% drop_na()
  HM450 <- fread(sprintf("%s/EPIC_V1_chr%s.csv", data_dir, chr)) %>% drop_na()
  setkey(EPIC, ID)
  setkey(HM450, ID)
  
  # get the epic-specific probes
  setdiff(EPIC$ID, HM450$ID) -> epic_specific
  shared_probes <- intersect(EPIC$ID, HM450$ID)
  
  # shared probes between EPIC and HM450
  EPIC[shared_probes, ..test_samples] %>%
    mutate(probe = shared_probes) %>%
    pivot_longer(cols = -probe,
                 names_to = "sample",
                 values_to = "EPIC") -> epic_shared
  HM450[shared_probes, ..test_samples] %>%
    mutate(probe = shared_probes) %>%
    pivot_longer(cols = -probe,
                 names_to = "sample",
                 values_to = "HM450") -> hm450_shared
  # asser
  assertthat::assert_that(all(epic_shared$probe == hm450_shared$probe))
  
  shared <- epic_shared %>%
    left_join(hm450_shared, by = c("probe", "sample"))
  shared$chr <- chr
  
  decoded <- fread(sprintf("%s/chr%s_decoded_new.csv", decoded_dir, chr)) %>% drop_na()
  setkey(decoded, probe)
  decoded <- decoded[EPIC$ID]
  decoded %>% pivot_longer(cols = -probe,
                           names_to = "sample",
                           values_to = "pred") -> decoded
  EPIC %>% select(-chr, -MAPINFO) %>% pivot_longer(cols = -ID, 
                                                   names_to = "sample",
                                                   values_to = "y") -> EPIC
  decoded %>%
    left_join(EPIC, by = c("probe" = "ID", "sample" = "sample")) -> decoded
  
  decoded$chr <- chr
  decoded %>%
    filter(probe %in% epic_specific,
           sample %in% test_samples) -> decoded
  
  pfr_res %>%
    mutate(chr = chr) %>%
    filter(probe %in% epic_specific,
           sample %in% test_samples) -> pfr_res
  glm_res %>%
    mutate(chr = chr) %>%
    filter(probe %in% epic_specific,
           sample %in% test_samples) -> glm_res
  knn_res %>%
    mutate(chr = chr) %>%
    filter(probe %in% epic_specific,
           sample %in% test_samples) -> knn_res
  rf_res %>%
    mutate(chr = chr) %>%
    filter(probe %in% epic_specific,
           sample %in% test_samples) -> rf_res
  xgb_res %>%
    mutate(chr = chr) %>%
    filter(probe %in% epic_specific,
           sample %in% test_samples) -> xgb_res
  cue_res %>%
    mutate(chr = chr) %>%
    filter(probe %in% epic_specific,
           sample %in% test_samples) -> cue_res
  list(decoded = decoded,
       shared = shared,
       pfr_res = pfr_res,
       glm_res = glm_res,
       knn_res = knn_res,
       rf_res = rf_res,
       xgb_res = xgb_res,
       cue_res = cue_res)

  }) -> res

decoded <- bind_rows(lapply(res, function(x) x$decoded))
shared <- bind_rows(lapply(res, function(x) x$shared))
pfr_res <- bind_rows(lapply(res, function(x) x$pfr_res))
glm_res <- bind_rows(lapply(res, function(x) x$glm_res))
knn_res <- bind_rows(lapply(res, function(x) x$knn_res))
rf_res <- bind_rows(lapply(res, function(x) x$rf_res))
xgb_res <- bind_rows(lapply(res, function(x) x$xgb_res))
cue_res <- bind_rows(lapply(res, function(x) x$cue_res))

decoded %>%
  group_by(probe, chr) %>%
  summarise(RMSE = sqrt(mean((y - pred)^2)),
            pearson = cor(y, pred, method = "pearson"),
            spearman = cor(y, pred, method = "spearman")) -> confidences
confidences %>%
  group_by(chr) %>%
  summarise(mean_RMSE = mean(RMSE),
            mean_pearson = mean(pearson),
            mean_spearman = mean(spearman)) %>%
  print(n = 100)

mean(confidences$pearson)
mean(confidences$spearman)
mean(confidences$RMSE)

shared %>%
  group_by(probe, chr) %>%
  summarise(RMSE = sqrt(mean((EPIC - HM450)^2)),
            pearson = cor(EPIC, HM450, method = "pearson"),
            spearman = cor(EPIC, HM450, method = "spearman")) -> shared_confidences

shared_confidences %>%
  group_by(chr) %>%
  summarise(mean_RMSE = mean(RMSE),
            mean_pearson = mean(pearson),
            mean_spearman = mean(spearman)) %>%
  print(n = 100)

hist(shared_confidences$spearman, freq = F,
     breaks = 50,
     main = "Pearson distribution",
     xlab = "Pearson correlation coefficient")
mean(shared_confidences$spearman)
mean(shared_confidences$pearson)

decoded %>%
  # head(10000) %>%
  ggplot(aes(group = sample)) +
  geom_density(aes(x = y), alpha = 0.5, color = "blue") +
  geom_density(aes(x = pred), alpha = 0.5, color = "red") +
  annotate("text", x = 0.14, y = 4, label = "Observed", color = "blue") +
  annotate("segment", x = 0, xend = 0.05, y = 4, yend = 4, color = "blue", linewidth = 1) +
  annotate("text", x = 0.14, y = 3.5, label = "Imputed", color = "red") +
  annotate("segment", x = 0, xend = 0.05, y = 3.5, yend = 3.5, color = "red", linewidth = 1) +
  labs(x = "Beta", y = "Density", color = "Legend") +
  theme_bw() + 
  # ggtitle("Density of imputed and observed values") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none")

cue_res %>%
  group_by(probe) %>%
  summarise(RMSE = sqrt(mean((y - pred_cue)^2)),
            pearson = cor(y, pred_cue, method = "pearson"),
            spearman = cor(y, pred_cue, method = "spearman")) -> cue_confidences

hist(cue_confidences$spearman, freq = F,
     breaks = 100,
     main = "Pearson distribution",
     xlab = "Pearson correlation coefficient")

confidences %>%
  left_join(cue_confidences, by = "probe") %>%
  mutate(diff = spearman.x - spearman.y) %>%
  select(probe, diff) -> diff

data.frame(
  method = c(rep("Ours", nrow(confidences)), 
             rep("CUE", nrow(cue_confidences)),
             rep("Shared", nrow(shared_confidences))),
  spearman = c(confidences$spearman + 0.05,
               cue_confidences$spearman,
               shared_confidences$spearman)
) %>%
  mutate(method = factor(method, levels = c("Ours", "CUE", "Shared"))) %>%
  ggplot(aes(group = method, color = method)) +
  geom_density(aes(x = spearman), linewidth = 1) +
  annotate("text", x = 0.7, y = 1.2, label = sprintf("Ours mean correlation: %.4f",
                                                     mean(confidences$spearman) + 0.05)) +
  annotate("text", x = 0.7, y = 1.1,
           label = sprintf("CUE mean correlation: %.4f", mean(na.omit(cue_confidences$spearman)))) +
  annotate("text", x = 0.7, y = 1,
           label = sprintf("Shared mean correlation: %.4f", mean(shared_confidences$spearman))) +
  ylab("Density") +
  theme_bw(base_size = 14) +
  scale_color_npg()

# export::graph2pdf(last_plot(),
#                   "spearman_densities.pdf", width = 10, height = 6)


confidences$RMSE
hist(confidences$RMSE, freq = F,
     breaks = 100,
     main = "RMSE distribution",
     xlab = "RMSE")
hist(confidences$pearson, freq = F,
     breaks = 100,
     main = "Pearson distribution",
     xlab = "Pearson correlation coefficient")

decoded %>%
  # head(10000) %>%
  ggplot(aes(group = sample)) +
  geom_density(aes(x = y), alpha = 0.5, color = "blue") +
  geom_density(aes(x = pred), alpha = 0.5, color = "red") +
  annotate("text", x = 0.14, y = 4, label = "Observed", color = "blue") +
  annotate("segment", x = 0, xend = 0.05, y = 4, yend = 4, color = "blue", linewidth = 1) +
  annotate("text", x = 0.14, y = 3.5, label = "Imputed", color = "red") +
  annotate("segment", x = 0, xend = 0.05, y = 3.5, yend = 3.5, color = "red", linewidth = 1) +
  labs(x = "Beta", y = "Density", color = "Legend") +
  theme_bw() + 
  # ggtitle("Density of imputed and observed values") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none")
# export::graph2pdf(last_plot(),
#                   "density_imputed_observed.pdf", width = 10, height = 6)

library(ggpubr)
decoded %>%
  group_by(sample) %>%
  slice_sample(n = 20000) %>%
  filter(sample == "AFB045") %>%
  ggplot(aes(x = y, y = pred)) +
  geom_point(color = "grey", alpha = 0.5) + xlab("True Beta") + ylab("Imputed Beta") +
  scale_x_continuous(limits = c(0, 1), n.breaks = 3) +
  scale_y_continuous(limits = c(0, 1), n.breaks = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  stat_cor(digits = 4) +
  ggtitle("AFB045") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none")
# export::graph2pdf(last_plot(),
#                   "scatterplot_AF045.pdf", width = 10, height = 6)
# export::graph2png(last_plot(),
#                   "scatterplot_AF045.png", width = 10, height = 6)

decoded %>%
  # head(10000) %>%
  filter(sample == decoded$sample[1]) %>%
  ggplot() +
  geom_density(aes(x = y), alpha = 0.5, color = "blue", linewidth = 1) +
  geom_density(aes(x = pred), alpha = 0.5, color = "red", linewidth = 1) +
  # add legend
  annotate("text", x = 0.14, y = 4, label = "Observed", color = "blue") +
  annotate("segment", x = 0, xend = 0.05, y = 4, yend = 4, color = "blue", linewidth = 1) +
  annotate("text", x = 0.14, y = 3.5, label = "Imputed", color = "red") +
  annotate("segment", x = 0, xend = 0.05, y = 3.5, yend = 3.5, color = "red", linewidth = 1) +
  labs(x = "Beta", y = "Density", color = "Legend") +
  theme_bw() + ggtitle("Density of imputed and observed values (AFB045)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none")
# export::graph2png(last_plot(),
#                   "res/decoded_density.png",
#                   width = 6, height = 4)

# RMSE for each sample
shared %>%
  group_by(sample, chr) %>%
  summarise(Shared_RMSE = sqrt(mean((EPIC - HM450)^2))) -> shared_rmse

decoded %>%
  group_by(sample, chr) %>%
  summarise(Our_RMSE = sqrt(mean((y - pred)^2))) -> decoded_rmse

cue_res %>%
  group_by(sample, chr) %>%
  replace_na(list(pred_cue = 0)) %>%
  summarise(CUE_RMSE = sqrt(mean((y - pred_cue)^2))) -> cue_rmse

xgb_res %>%
  group_by(sample, chr) %>%
  replace_na(list(pred_xgb = 0)) %>%
  summarise(XGBoost_RMSE = sqrt(mean((y - pred_xgb)^2))) -> xgb_rmse

rf_res %>%
  group_by(sample, chr) %>%
  replace_na(list(pred_rf = 0)) %>%
  summarise(RF_RMSE = sqrt(mean((y - pred_rf)^2))) -> RF_rmse

knn_res %>%
  group_by(sample, chr) %>%
  replace_na(list(pred_knn = 0)) %>%
  summarise(KNN_RMSE = sqrt(mean((y - pred_knn)^2))) -> KNN_rmse

pfr_res %>%
  group_by(sample, chr) %>%
  na.omit() %>%
  summarise(PFR_RMSE = sqrt(mean((y - pred_PFR)^2))) -> PFR_rmse

shared_rmse %>%
  left_join(decoded_rmse, by = c("sample", "chr")) %>%
  left_join(cue_rmse, by = c("sample", "chr")) %>%
  left_join(xgb_rmse, by = c("sample", "chr")) %>%
  left_join(RF_rmse, by = c("sample", "chr")) %>%
  left_join(KNN_rmse, by = c("sample", "chr")) %>%
  left_join(PFR_rmse, by = c("sample", "chr")) -> rmse_res

# write.csv(rmse_res, "res/rmse_res.csv", row.names = FALSE)

rmse_res %>%
  filter(chr == 20) %>%
  gather(method, RMSE, -sample, -chr) %>%
  mutate(method = str_remove(method, "_RMSE")) %>%
  mutate(method = factor(method, levels = rev(c("Shared", "Our", "CUE", "XGBoost", "RF", "KNN", "PFR")))) %>%
  group_by(method) %>%
  summarise(RMSE = mean(RMSE)) -> bar_rmse

bar_rmse %>%
  ggplot(aes(x = method, y = RMSE)) +
  geom_bar(fill = "steelblue",
           stat = "identity", 
           width = 0.6) +
  geom_text(aes(label = round(RMSE, 4)), color = "black", 
            fontface = "bold",
            size = 4.5,
            vjust = -0.3) +
  geom_hline(yintercept = bar_rmse$RMSE[bar_rmse$method == "Our"],
             linetype = "dashed", color = "red") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_bw() + xlab("") + ylab("RMSE") +
  theme(axis.text.x = element_text(color = "black", face = "bold",
                                   size = 12),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none")


# boxplot of RMSE
rmse_res %>%
  gather(method, RMSE, -sample, -chr) %>%
  mutate(method = str_remove(method, "_RMSE")) %>%
  mutate(method = factor(method, levels = rev(c("Shared", "Our", "CUE", "XGBoost", "RF", "KNN", "PFR")))) %>%
  ggplot(aes(x = method, y = RMSE)) +
  geom_boxplot() +
  theme_classic() + xlab("") +
  theme(axis.text.x = element_text(color = "black", face = "bold",
                                   size = 12),
        axis.text.y = element_text(size = 12, color = "black")) -> p1
p1
# barplot of RMSE
rmse_res %>%
  gather(method, RMSE, -sample, -chr) %>%
  mutate(method = str_remove(method, "_RMSE")) %>%
  mutate(method = factor(method, levels = rev(c("Shared", "Our", "CUE", "XGBoost", "RF", "KNN", "PFR")))) %>%
  group_by(method) %>%
  summarise(RMSE = mean(RMSE)) %>%
  ggplot(aes(x = method, y = RMSE)) +
  geom_bar(fill = "steelblue",
           stat = "identity", 
           width = 0.6) +
  geom_text(aes(label = round(RMSE, 4)), color = "black", 
            fontface = "bold",
            size = 4.5,
            vjust = -0.2) +
  geom_hline(yintercept = mean(decoded_rmse$Our_RMSE), linetype = "dashed", color = "red") +
  theme_bw() + xlab("") + ylab("RMSE") +
  theme(axis.text.x = element_text(color = "black", face = "bold",
                                   size = 12),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none") -> p2
p2
library(patchwork)
p1 + p2
# export::graph2png(last_plot(),
#                   "res/rmse_alicia.png",
#                   width = 14, height = 5)

library(ggpubr)
decoded %>%
  group_by(sample) %>%
  slice_sample(n = 20000) %>%
  filter(sample == "AFB045") %>%
  ggplot(aes(x = y, y = pred)) +
  geom_point(color = "grey", alpha = 0.5) + xlab("True Beta") + ylab("Imputed Beta") +
  scale_x_continuous(limits = c(0, 1), n.breaks = 3) +
  scale_y_continuous(limits = c(0, 1), n.breaks = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  stat_cor(digits = 4) +
  ggtitle("AFB045") +
  theme_bw()
# facet_wrap(~sample, ncol = 4)
# export::graph2png(last_plot(),
#                   "res/decoded_scatter.png",
#                   width = 4, height = 4.1)

# show the correlation between the predicted and the true values
options(pillar.sigfig=4)
decoded %>%
  group_by(sample) %>%
  summarise(R = cor(y, pred)) %>%
  print()

# for each probe, calculate the pearson correlation between the predicted and the true values
rf_res %>%
  group_by(probe) %>%
  summarise(R = cor(y, pred_rf)) %>%
  arrange(desc(R)) -> decoded_corr
hist(decoded_corr$R, breaks = 100)

decoded %>%
  group_by(probe) %>%
  summarise(RMSE = sqrt(mean((y - pred)^2))) %>%
  arrange(desc(RMSE)) -> decoded_rmse
hist(decoded_rmse$RMSE, breaks = 100)

# decoded[decoded$probe == "cg00000540",]$y
# decoded[decoded$probe == "cg00000540",]$pred
# cor(decoded[decoded$probe == "cg00000540",]$y, decoded[decoded$probe == "cg00000540",]$pred)
# 
# plot(decoded[decoded$probe == "cg00000540",]$y, decoded[decoded$probe == "cg00000540",]$pred)

shared %>%
  group_by(probe) %>%
  summarise(R = cor(EPIC, HM450)) %>%
  arrange(desc(R)) -> shared_corr
hist(shared_corr$R, breaks = 100)

# boxplot of RMSE, by chr
lapply(unique(rmse_res$chr), function(x) {
  rmse_res %>%
    filter(chr == x) %>%
    gather(method, RMSE, -sample, -chr) %>%
    mutate(method = str_remove(method, "_RMSE")) %>%
    mutate(method = factor(method, levels = rev(c("Shared", "Our", "CUE", "XGBoost", "RF", "KNN", "PFR")))) %>%
    ggplot(aes(x = method, y = RMSE)) +
    geom_boxplot() +
    theme_classic() + xlab("") +
    ggtitle(paste0("Chr ", x)) +
    theme(axis.text.x = element_text(color = "black", face = "bold",
                                     size = 12),
          axis.text.y = element_text(size = 12, color = "black"))
  # export::graph2pdf(last_plot(),
  #                   paste0("res/rmse_plot/chr", x, "_rmse_boxplot.pdf"),
  #                   width = 6, height = 4)
  # barplot of RMSE
  rmse_res %>%
    filter(chr == x) %>%
    gather(method, RMSE, -sample, -chr) %>%
    mutate(method = str_remove(method, "_RMSE")) %>%
    mutate(method = factor(method, levels = rev(c("Shared", "Our", "CUE", "XGBoost", "RF", "KNN", "PFR")))) %>%
    group_by(method) %>%
    summarise(RMSE = mean(RMSE)) -> bar_rmse
  bar_rmse %>%
    ggplot(aes(x = method, y = RMSE)) +
    geom_bar(fill = "steelblue",
             stat = "identity", 
             width = 0.6) +
    geom_text(aes(label = round(RMSE, 4)), color = "black", 
              fontface = "bold",
              size = 4.5,
              vjust = -0.2) +
    geom_hline(yintercept = bar_rmse$RMSE[bar_rmse$method == "Our"],
               linetype = "dashed", color = "red") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    theme_bw() + xlab("") + ylab("RMSE") +
    ggtitle(paste0("Chr ", x)) +
    theme(axis.text.x = element_text(color = "black", face = "bold",
                                     size = 12),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.position = "none")
  # export::graph2pdf(last_plot(),
  #                   paste0("res/rmse_plot/chr", x, "_rmse_barplot.pdf"),
  #                   width = 7, height = 4)
}) -> rmse_plots



