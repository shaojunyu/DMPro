library(data.table)

data_dir <- "tmp/AnkeHuels_405K_EPIC/"
baseline_dir <- "res/AnkeHuels_405K_EPIC_baseline/"
decoded_dir <- "res/AnkeHuels_405K_EPIC_decoded/"

train_samples <- read.csv(paste0(data_dir, "/train_samples.txt"), header = F)$V1
test_samples <- read.csv(paste0(data_dir, "/test_samples.txt"), header = F)$V1
val_samples <- read.csv(paste0(data_dir, "/val_samples.txt"), header = F)$V1

HM450_chr20 <-
  fread("tmp/AnkeHuels_405K_EPIC/HM450_chr20.csv")

HM450_chr20 %>%
  select(chr, ID, AFB180) -> test_data_for_shiny
write.csv(test_data_for_shiny, "src/shiny_app/test_data_for_shiny.csv", row.names = F)
