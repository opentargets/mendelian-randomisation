# Extract MR results

## Load library and authenticate

source("loadLibraries_and_authenticate.r")

## list of names of protein datasets

n <- c("SUHRE_2017", "FOLKERSEN_2017", "OLLI_2017", "PIETZNER_2020", "HILLARY_2019", "DRAISMA_2015", "KETTUNEN_2016", "SHIN_2014")

## location of outputs

url <- "gs://genetics-portal-ukbb-mr-eur/forward-pan-gsmr-outputs/"

system("mkdir mr_results && mkdir mr_results/full && mkdir mr_results/filtered")

parallel::mclapply(seq(n), function (x) {
  
  dfn <- strsplit(n, "_")[[x]][1]
  
  gcs_get_object(paste0(url, dfn, "/results.gsmr"), saveToDisk = paste0(dfn, "_gsmr"))
  
  df <- fread(paste0(dfn, "_gsmr"), stringsAsFactors = F)
  
  df[df=="NaN"] <- NA_character_
  
  df2 <- na.omit(df)
  
  df2$Data <- n[x]
  
  df2$p_fdr <- p.adjust(df2$p, method = "BH")
  
  df2$p <- as.numeric(df2$p)
  
  # df2 <- df2[df2$p < 0.05,] # uncomment this and next line to extract filtered dataset
  
  # write.csv(df2, paste0("mr_results/filtered/", n[x], "_mr_filtered.csv"), row.names = F, quote = F)
  
  write.csv(df2, paste0("mr_results/full/", n[x], "_mr_full.csv"), row.names = F, quote = F)
  
  system(paste0("rm ", dfn, "_gsmr"))
  
}, mc.cores = parallel::detectCores())