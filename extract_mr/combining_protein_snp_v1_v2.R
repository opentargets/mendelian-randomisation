# combining full and full_v2 for proteins

library(dplyr)

prot <- paste(c("SUN", "PIETZNER", "SCALLOP", "SUHRE", "FOLKERSEN", "HILLARY", "OLLI"), collapse = "|")

v1 <- list.files("~/mr_results_snps/full", full.names = T, pattern = prot) 

v2 <- list.files("~/mr_results_snps/full_v2", full.names = T, pattern = prot)

f <- function (x) {
  
  data.table::fread(x, stringsAsFactors = F)
  
}

df <- rbindlist(lapply(v1, f))

df2 <- rbindlist(lapply(v2, f))

df3 <- rbindlist(list(df, df2))

saveRDS(df3, "~/mr_results_snps/full_protein_MR_snp_v1_v2_feb2022.rds")
