# Extract MR results per SNP

# GSMR per-SNP results are stored in results.eff_plot.gz that takes a lot of time to extract.

# This script 1) extracts the full data from results.eff_plot.gz using a function from gsmr_plot_v2, and 2) stores the full SNP data as an RDS object so it can be loaded easily and queried

## Generating RDS objects

prot <- c("SUHRE", "FOLKERSEN", "HILLARY", "OLLI", "DRAISMA", "SHIN")

dirs <- gcs_list_objects(bucket = "genetics-portal-ukbb-mr-eur", prefix = paste0("forward-pan-gsmr-outputs/", prot))

source("gsmr_plot_v2.r") # hacked the gsmr_snp_effect function to show MR results with less than 3 instruments

parallel::mclapply(seq(prot), function (x) {

  system(paste0("mkdir ", prot[x], "_gsmr_data"))

  gcs_get_object(paste0("gs://genetics-portal-ukbb-mr-eur/", "forward-pan-gsmr-outputs/", prot[x], "/results.eff_plot.gz"), saveToDisk = paste0(prot[x], "_gsmr_data/", prot[x], "_results_gsmr_snp.gz"))

  mrdf <- read_gsmr_data(paste0(prot[x], "_gsmr_data/", prot[x], "_results_gsmr_snp.gz")) # this takes a bit of time depending on the size of the dataset

  saveRDS(mrdf, paste0(prot[x], "_gsmr_data/", prot[x], "_results_gsmr_snp.rds"))

}, mc.cores = parallel::detectCores())

## Using the full MR results on RDS objects to extracts SNPs used for every MR associations

n <- c("SUHRE", "FOLKERSEN", "HILLARY", "OLLI", "DRAISMA", "SHIN")


parallel::mclapply(seq(n), function (x) {
  
  # load per-SNP rds objects which are lists
  
  url <- "gs://genetics-portal-ukbb-mr-eur/protein-gsmr-data/"
  
  n2 <- strsplit(n, "_")[[x]][1]
  
  dfn <- paste0(n2, "_gsmr_data/", n2, "_results_gsmr_snp.rds")
  
  gcs_get_object(paste0(url, dfn), saveToDisk = n2)
  
  dfsnp <- readRDS(n2)
  
  # load full MR results
  
  dfmr <- data.table::fread(paste0("mr_results/full/", n[x], "_mr_full.csv"), stringsAsFactors = F)
  
  genes <- dfmr$Exposure
  
  out <- dfmr$Outcome
  
  bxy <- dfmr$bxy
  
  bxy_se <- dfmr$se
  
  bxy_pval <- dfmr$p
  
  nsnps <- dfmr$nsnp
  
  Data <- dfmr$Data
  
  # Use full MR results to extract per-SNP associations for each dataset
  
  dflst <- lapply(seq_along(genes), function (y) {
    
    # print(genes[y])
    
    df3 <- gsmr_snp_effect(gsmr_data = dfsnp, expo_str = genes[y], outcome_str = out[y])
    
    df4 <- data.frame(df3, stringsAsFactors = F)
    
    df4$bxy <- bxy[y]
    
    df4$bxy_se <- bxy_se[y]
    
    df4$bxy_pval <- bxy_pval[y]
    
    df4$exposure <- genes[y]
    
    df4$outcome <- out[y]
    
    df4$nsnp <- nsnps[y]
    
    df4$Data <- Data[y]
    
    df4
    
  })
  
  df3 <- rbindlist(dflst)
  
  write.csv(df3, paste0("mr_results_snps/full/", n[x], "_mr_full_snps.csv"), row.names = F, quote = F)
  
  system(paste0("rm ", n2))
  
}, mc.cores = parallel::detectCores())

