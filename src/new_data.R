library(minfi)
library(minfiData)

HM450K_sample <- openxlsx::read.xlsx(
  "data/Steve_idats/EPICv1v2_sampleList.xlsx",
  sheet = 2)
EPCI_sample <- openxlsx::read.xlsx(
  "data/Steve_idats/EPICv1v2_sampleList.xlsx",
  sheet = 1
)

baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)
