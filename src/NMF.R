library(data.table)

V1 <- fread("tmp/EPIC_V1_V2/EPIC_V1_chr20.csv")
V2 <- fread("tmp/EPIC_V1_V2/EPIC_V2_chr20.csv")

# remove duplicated probes
V1 <- V1[!duplicated(V1$ID),]
V2 <- V2[!duplicated(V2$ID),]


shared_probes <- intersect(V1$ID, V2$ID)
V2_specific_probes <- setdiff(V2$ID, V1$ID)

length(unique(V2$ID)) - length(unique(V1$ID))

V1_shared <- V1[V1$ID %in% shared_probes,]
V2_shared <- V2[V2$ID %in% shared_probes,]
plot(V1_shared$`12123-C`, V2_shared$`12123-C`)

plot(c(V1_shared[1,4:47]), c(V2_shared[1,4:47]))
cor(as.numeric(V1_shared[1,4:47]), as.numeric(V2_shared[1,4:47]))

lapply(1:nrow(V1_shared), function(i) {
  cor(as.numeric(V1_shared[i,4:47]), as.numeric(V2_shared[i,4:47]))
}) -> cor_values
hist(unlist(cor_values))
mean(unlist(cor_values))

hist(V1_shared$`12123-C`)
hist(V2_shared$`12123-C`)

V2_test <- V2[V2$ID %in% c(shared_probes,V2_specific_probes),]







