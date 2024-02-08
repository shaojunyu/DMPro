# This script is used to preprocess the raw data from 450k and EPIC arrays.
# For AliciaKSmith_450K_EPIC data

library(tidyverse)
library(minfi)

EPIC_450K_sample <- openxlsx::read.xlsx(
  "data/AliciaKSmith_450K_EPIC/Steve_idats/EPICv1v2_sampleList.xlsx",
  sheet = 2)

HM450K_basenames <- 
  sprintf("data/AliciaKSmith_450K_EPIC/Steve_idats/450k_realSamples/%s",
          EPIC_450K_sample$ID_450k_blood)

HM450_RGSet <- read.metharray(HM450K_basenames, verbose = T)
# preprocessing
HM450_Mset <- preprocessIllumina(HM450_RGSet)
HM450_beta <- getBeta(HM450_Mset)

list.files("data/AliciaKSmith_450K_EPIC/GSE132203", full.names = T) %>%
  lapply(function(x){
    basename(x) %>% 
      str_split("_", simplify = T) %>%
      .[,2:3] %>%
      paste(collapse = "_") -> tmp
    if (tmp %in% EPIC_450K_sample$EPICv1) {
      x
    }else{
      file.remove(x)
    }
  }) %>% unlist() -> valid_EPIC_basenames

valid_EPIC_basenames %>%
  str_remove("_Red.idat.gz") %>%
  str_remove("_Grn.idat.gz") %>%
  unique() -> valid_EPIC_basenames

EPIC_RGSet <- read.metharray(valid_EPIC_basenames, verbose = T, force = T)
EPIC_MSet <- preprocessIllumina(EPIC_RGSet)
EPIC_beta <- getBeta(EPIC_MSet)

colnames(EPIC_beta) %>%
  str_remove("GSM\\d+_") %in% EPIC_450K_sample$EPICv1

colnames(EPIC_beta) %>%
  str_remove("GSM\\d+_") -> colnames(EPIC_beta)

dim(HM450_beta)
dim(EPIC_beta)

head(HM450_beta)
head(EPIC_beta)

colnames(HM450_beta) %in% EPIC_450K_sample$ID_450k_blood
colnames(EPIC_beta) %in% EPIC_450K_sample$EPICv1

HM450_beta[,EPIC_450K_sample$ID_450k_blood] -> HM450_beta
EPIC_beta[,EPIC_450K_sample$EPICv1] -> EPIC_beta
colnames(HM450_beta) <- EPIC_450K_sample$Sample_Name
colnames(EPIC_beta) <- EPIC_450K_sample$Sample_Name

dim(na.omit(HM450_beta))
dim(na.omit(EPIC_beta))

# clean data
na.omit(HM450_beta) -> HM450_beta
na.omit(EPIC_beta) -> EPIC_beta

# read annotation
HM450K_annotation <- getAnnotation(HM450_Mset)
EPIC_annotation <- getAnnotation(EPIC_MSet)
rownames(HM450_beta)
rownames(EPIC_beta)

