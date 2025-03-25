# merge Anke's data and Alicia's data

Anke_450K <- data.table::fread("data/AnkeHuels_405K_EPIC/HM450_beta.csv")
Anke_EPIC <- data.table::fread("data/AnkeHuels_405K_EPIC/EPIC_beta.csv")

Alicia_450K <- data.table::fread("data/AliciaKSmith_450K_EPIC/HM450_beta.csv", header = T)
Alicia_EPIC <- data.table::fread("data/AliciaKSmith_450K_EPIC/EPIC_beta.csv", header = T)
