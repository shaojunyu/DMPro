library(data.table)
library(tidyverse)

# Evaluate cross-dataset (transfer) imputation across all 22 chromosomes.
# Compares decoded predictions against true EPIC values on EPIC-specific probes
# (probes on EPIC but not 450K) — the ones actually being imputed.

data_dir    <- "tmp/AliciaKSmith_450K_EPIC"
decoded_dir <- "res/AliciaKSmith_450K_EPIC_decoded_AH"
out_dir     <- "res/AliciaKSmith_450K_EPIC_plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

per_sample_chr <- list()

for (chr in 1:22) {
  decoded_file <- sprintf("%s/chr%s_decoded.csv", decoded_dir, chr)
  if (!file.exists(decoded_file)) {
    message(sprintf("chr%s: no decoded file, skipping", chr))
    next
  }

  res <- tryCatch({
    EPIC  <- fread(sprintf("%s/EPIC_chr%s.csv",  data_dir, chr)) %>% drop_na()
    HM450 <- fread(sprintf("%s/HM450_chr%s.csv", data_dir, chr)) %>% drop_na()
    setkey(EPIC, ID)
    setkey(HM450, ID)

    # EPIC-specific probes = the imputation target
    epic_specific <- setdiff(EPIC$ID, HM450$ID)

    decoded <- fread(decoded_file) %>% drop_na()
    setkey(decoded, probe)
    decoded <- decoded[EPIC$ID]
    decoded <- decoded %>%
      pivot_longer(cols = -probe, names_to = "sample", values_to = "pred")

    sample_names <- unique(decoded$sample)
    EPIC_long <- EPIC %>%
      select(ID, all_of(sample_names)) %>%
      pivot_longer(cols = -ID, names_to = "sample", values_to = "y")

    decoded <- decoded %>%
      left_join(EPIC_long, by = c("probe" = "ID", "sample" = "sample")) %>%
      filter(probe %in% epic_specific)
    decoded$chr <- chr

    decoded %>%
      drop_na(y, pred) %>%
      group_by(sample, chr) %>%
      summarise(Our_RMSE = sqrt(mean((y - pred)^2)),
                n_probes = n(), .groups = "drop")
  }, error = function(e) {
    message(sprintf("chr%s: error — %s", chr, conditionMessage(e)))
    NULL
  })

  if (!is.null(res)) {
    per_sample_chr[[as.character(chr)]] <- res
    message(sprintf("chr%s: done (%d samples, mean RMSE %.4f)",
                    chr, nrow(res), mean(res$Our_RMSE)))
  }
}

per_sample_chr <- bind_rows(per_sample_chr)

if (nrow(per_sample_chr) == 0) {
  stop("No results produced — check that decoded files exist in ", decoded_dir)
}

# per-chromosome summary (simple mean of per-sample RMSE)
per_chr <- per_sample_chr %>%
  group_by(chr) %>%
  summarise(mean_RMSE = mean(Our_RMSE),
            n_samples = n(),
            n_probes  = sum(n_probes), .groups = "drop") %>%
  arrange(chr)

# overall summary (simple mean of per-sample RMSE across all chr & samples)
overall <- per_sample_chr %>%
  summarise(mean_RMSE = mean(Our_RMSE),
            n_probes  = sum(n_probes))

fwrite(per_sample_chr, file.path(out_dir, "transfer_rmse_per_sample_chr.csv"))
fwrite(per_chr,        file.path(out_dir, "transfer_rmse_per_chr.csv"))

cat("\n=== Per-chromosome RMSE (mean of per-sample RMSE) ===\n")
print(per_chr, n = 22)
cat(sprintf("\n=== Overall transfer RMSE (EPIC-specific probes): %.4f  (%d probe-values) ===\n",
            overall$mean_RMSE, overall$n_probes))
