library(data.table)
library(tidyverse)

chr <- 20
data_dir <- "tmp/AliciaKSmith_450K_EPIC/"
decoded_dir <- "res/AliciaKSmith_450K_EPIC_decoded_AH/"

# read data
EPIC <- fread(sprintf("%s/EPIC_chr%s.csv", data_dir, chr)) %>% drop_na()
HM450 <- fread(sprintf("%s/HM450_chr%s.csv", data_dir, chr)) %>% drop_na()
setkey(EPIC, ID)
setkey(HM450, ID)

# get the epic-specific probes
setdiff(EPIC$ID, HM450$ID) -> epic_specific
shared_probes <- intersect(EPIC$ID, HM450$ID)

decoded <- fread(sprintf("%s/chr%s_decoded.csv", decoded_dir, chr)) %>% drop_na()
setkey(decoded, probe)
decoded <- decoded[EPIC$ID]
decoded %>% pivot_longer(cols = -probe,
                         names_to = "sample",
                         values_to = "pred") -> decoded
sample_names <- unique(decoded$sample)
EPIC %>%
  select(ID, chr, MAPINFO, all_of(sample_names)) %>%
  select(-chr, -MAPINFO) %>% pivot_longer(cols = -ID, 
                                          names_to = "sample",
                                          values_to = "y") -> EPIC
decoded %>%
  left_join(EPIC, by = c("probe" = "ID", "sample" = "sample")) -> decoded
decoded$chr <- chr

decoded %>%
  filter(probe %in% epic_specific) -> decoded

decoded %>%
  group_by(sample, chr) %>%
  summarise(Our_RMSE = sqrt(mean((y - pred)^2)))
