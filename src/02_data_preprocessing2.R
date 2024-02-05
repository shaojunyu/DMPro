# This script is used to preprocess the raw data from 450k and EPIC arrays.
# For AliciaKSmith_450K_EPIC data

library(minfi)
library(minfiData)

EPIC_450K_sample <- openxlsx::read.xlsx(
  "data/AliciaKSmith_450K_EPIC/Steve_idats/EPICv1v2_sampleList.xlsx",
  sheet = 2)
