# The per-SNP parts depend on the rds files. So extract the effect.gz files and create rds files (extracting_per_snp.rds) first.

library(googleCloudStorageR)
gcs_auth(json_file = "open-targets-genetics-fc5b6cda58e5.json")

library(data.table)

ver <- "v2"

# MR TABLES (FULL)====

# ALL EXCCEPT SUN AND NIGHTINGALE

n <- c("SCALLOP_2020", "SUHRE_2017", "FOLKERSEN_2017", "OLLI_2017", "PIETZNER_2020", "HILLARY_2019", "DRAISMA_2015", "KETTUNEN_2016", "SHIN_2014")

system(paste0("mkdir ", "mr_results/full_", ver))

parallel::mclapply(seq(n), function (x) {
  
  url <- paste0("gs://genetics-portal-ukbb-mr-eur/forward-pan-gsmr-outputs-nsnp1_", ver, "/")
  
  dfn <- strsplit(n, "_")[[x]][1]
  
  gcs_get_object(paste0(url, dfn, "/results.gsmr"), saveToDisk = paste0(dfn, "_gsmr"))
  
  df <- fread(paste0(dfn, "_gsmr"), stringsAsFactors = F)
  
  df[df=="NaN"] <- NA_character_
  
  df2 <- na.omit(df)
  
  df2$Data <- n[x]
  
  write.csv(df2, paste0("mr_results/full_", ver, "/", n[x], "_mr_full.csv"), row.names = F, quote = F)
  
  system(paste0("rm ", dfn, "_gsmr"))
  
}, mc.cores = parallel::detectCores())

# SUN and NIGHTINGALE

n <- c("SUN_2018", "NIGHTINGALE_2021")

ver <- "v2"

a <- seq(1,11,1)

dir <- sprintf("dir_%03d", a)

lapply(seq(n), function (y) { # extracting the files from the directory
  
  dfn <- strsplit(n, "_")[[y]][1]
  
  system(paste0("mkdir ", "mr_results/", dfn ,"_", ver, "/"))
  
  parallel::mclapply(seq(dir), function (x) {
    
    d <- dir[x]
    
    url <- paste0("gs://genetics-portal-ukbb-mr-eur/forward-pan-gsmr-outputs-nsnp1_", ver, "/", dfn, "/")
                  
                  gcs_get_object(paste0(url, d, "/results.gsmr"), saveToDisk = paste0(dfn, "_", d, "_gsmr"))
                  
                  df <- fread(paste0(dfn, "_", d, "_gsmr"), stringsAsFactors = F)
                  
                  df[df=="NaN"] <- NA_character_
                  
                  df2 <- na.omit(df)
                  
                  df2$Data <- n[y]
                  
                  write.csv(df2, paste0("mr_results/", dfn ,"_", ver, "/", d, "_mr_full.csv"), row.names = F, quote = F)
                  
                  system(paste0("rm ", dfn, "_", d, "_gsmr"))
                  
  }, mc.cores = parallel::detectCores())
    
  
})


lapply(seq(n), function (y) { # combining the extracted files to a single csv file
  
  dfn <- strsplit(n, "_")[[y]][1]
  
  print(dfn)
  
  p <- paste0("~/mr_results/", dfn, "_", ver)
  
  # df <- do.call(rbind, lapply(list.files(path = p), read.csv))
  
  f <- function (x) {
    
    data.table::fread(x, stringsAsFactors = F)
    
  }
  
  fn <- list.files(path = p, full.names = T)
  
  print(fn)
  
  df <- rbindlist(lapply(fn, f))
  
  print(names(df))
  
  p <- paste0("~/mr_results/full_", ver, "/", n[y], "_mr_full.csv")
  
  print(p)
  
  write.csv(df, p, row.names = F, quote = F)
  
})



# MR TABLES PER SNP (FULL)====

# ALL EXCCEPT SUN AND NIGHTINGALE

n <- c("SCALLOP_2020", "SUHRE_2017", "FOLKERSEN_2017", "OLLI_2017", "PIETZNER_2020", "HILLARY_2019", "DRAISMA_2015", "KETTUNEN_2016", "SHIN_2014")

system(paste0("mkdir ", "mr_results_snps/full_", ver))

parallel::mclapply(seq(n), function (x) {
  
  # load per-SNP rds objects which are lists
  
  # url <- "gs://genetics-portal-ukbb-mr-eur/protein-gsmr-data/"
  
  n2 <- strsplit(n, "_")[[x]][1]
  
  dfn <- paste0(n2, "_gsmr_", ver, "_data/", n2, "_results_gsmr_snp.rds")
  
  # gcs_get_object(paste0(url, dfn), saveToDisk = n2)
  
  dfsnp <- readRDS(dfn)
  
  # load full MR results
  
  dfmr <- data.table::fread(paste0("~/mr_results/full_", ver, "/", n[x], "_mr_full.csv"), stringsAsFactors = F)
  
  genes <- dfmr$Exposure
  
  out <- dfmr$Outcome
  
  bxy <- dfmr$bxy
  
  bxy_se <- dfmr$se
  
  bxy_pval <- dfmr$p
  
  nsnps <- dfmr$nsnp
  
  Data <- dfmr$Data
  
  source("gsmr_plot.r")
  
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
  
  write.csv(df3, paste0("~/mr_results_snps/full_", ver, "/", n[x], "_mr_full_snps.csv"), row.names = F, quote = F)
  
  # system(paste0("rm ", n2))
  
}, mc.cores = parallel::detectCores())

# SUN AND NIGHTINGALE

n <- c("SUN_2018", "NIGHTINGALE_2021")

ver <- "v2"

a <- seq(1,11,1)

dir <- sprintf("dir_%03d", a)

source("gsmr_plot.r")

lapply(seq(n), function (z) {
  
  n2 <- strsplit(n, "_")[[z]][1]
  
  system(paste0("mkdir ", "mr_results_snps/", n2 ,"_", ver, "/"))
  
  parallel::mclapply(seq(dir), function (x) {
    
    require(data.table)
    
    # load per-SNP rds objects which are lists
    
    # url <- "gs://genetics-portal-ukbb-mr-eur/protein-gsmr-data/"
    
    # n2 <- strsplit(n, "_")[[x]][1]
    
    dfn <- paste0(n2, "_gsmr_", ver, "_data/", dir[x], "_results_gsmr_snp.rds")
    
    # gcs_get_object(paste0(url, dfn), saveToDisk = n2)
    
    dfsnp <- readRDS(dfn)
    
    # load full MR results
    
    dfmr <- data.table::fread(paste0("~/mr_results/", n2, "_", ver, "/", dir[x], "_mr_full.csv"), stringsAsFactors = F)
    
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
    
    df3 <- data.table::rbindlist(dflst)
    
    write.csv(df3, paste0("~/mr_results_snps/", n2, "_", ver, "/", dir[x], "_mr_full_snps.csv"), row.names = F, quote = F)
    
    # system(paste0("rm ", n2))
    
  }, mc.cores = parallel::detectCores())
  
  
})


lapply(seq(n), function (y) {
  
  dfn <- strsplit(n, "_")[[y]][1]
  
  print(dfn)
  
  p <- paste0("~/mr_results_snps/", dfn, "_", ver)
  
  f <- function (x) {
    
    data.table::fread(x, stringsAsFactors = F)
    
  }
  
  fn <- list.files(path = p, full.names = T)
  
  print(fn)
  
  df <- rbindlist(lapply(fn, f))
  
  p <- paste0("~/mr_results_snps/full_", ver, "/", n[y], "_mr_full.csv")
  
  print(p)
  
  write.csv(df, p, row.names = F, quote = F)
  
})

# Special run for Nightingale v1

## MR table

n <- "NIGHTINGALE_2021"

# ver <- "v2"

a <- seq(1,11,1)

dir <- sprintf("dir_%03d", a)

lapply(seq(n), function (y) { # extracting the files from the directory
  
  dfn <- strsplit(n, "_")[[y]][1]
  
  system(paste0("mkdir ", "mr_results/", dfn , "/"))
  
  parallel::mclapply(seq(dir), function (x) {
    
    d <- dir[x]
    
    url <- paste0("gs://genetics-portal-ukbb-mr-eur/forward-pan-gsmr-outputs-nsnp1/", dfn, "/")
    
    gcs_get_object(paste0(url, d, "/results.gsmr"), saveToDisk = paste0(dfn, "_", d, "_gsmr"))
    
    df <- fread(paste0(dfn, "_", d, "_gsmr"), stringsAsFactors = F)
    
    df[df=="NaN"] <- NA_character_
    
    df2 <- na.omit(df)
    
    df2$Data <- n[y]
    
    write.csv(df2, paste0("mr_results/", dfn, "/", d, "_mr_full.csv"), row.names = F, quote = F)
    
    system(paste0("rm ", dfn, "_", d, "_gsmr"))
    
  }, mc.cores = parallel::detectCores())
  
  
})


lapply(seq(n), function (y) { # combining the extracted files to a single csv file
  
  dfn <- strsplit(n, "_")[[y]][1]
  
  print(dfn)
  
  p <- paste0("~/mr_results/", dfn)
  
  # df <- do.call(rbind, lapply(list.files(path = p), read.csv))
  
  f <- function (x) {
    
    data.table::fread(x, stringsAsFactors = F)
    
  }
  
  fn <- list.files(path = p, full.names = T)
  
  print(fn)
  
  df <- rbindlist(lapply(fn, f))
  
  print(names(df))
  
  p <- paste0("~/mr_results/full/", n[y], "_mr_full.csv")
  
  print(p)
  
  write.csv(df, p, row.names = F, quote = F)
  
})

## MR SNP table

source("gsmr_plot.r")

lapply(seq(n), function (z) {
  
  n2 <- strsplit(n, "_")[[z]][1]
  
  system(paste0("mkdir ", "mr_results_snps/", n2, "/"))
  
  parallel::mclapply(seq(dir), function (x) {
    
    require(data.table)
    
    # load per-SNP rds objects which are lists
    
    # url <- "gs://genetics-portal-ukbb-mr-eur/protein-gsmr-data/"
    
    # n2 <- strsplit(n, "_")[[x]][1]
    
    dfn <- paste0(n2, "_gsmr_data/", dir[x], "_results_gsmr_snp.rds")
    
    # gcs_get_object(paste0(url, dfn), saveToDisk = n2)
    
    dfsnp <- readRDS(dfn)
    
    # load full MR results
    
    dfmr <- data.table::fread(paste0("~/mr_results/", n2, "/", dir[x], "_mr_full.csv"), stringsAsFactors = F)
    
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
    
    df3 <- data.table::rbindlist(dflst)
    
    write.csv(df3, paste0("~/mr_results_snps/", n2, "/", dir[x], "_mr_full_snps.csv"), row.names = F, quote = F)
    
    # system(paste0("rm ", n2))
    
  }, mc.cores = parallel::detectCores())
  
  
})


lapply(seq(n), function (y) {
  
  dfn <- strsplit(n, "_")[[y]][1]
  
  print(dfn)
  
  p <- paste0("~/mr_results_snps/", dfn)
  
  f <- function (x) {
    
    data.table::fread(x, stringsAsFactors = F)
    
  }
  
  fn <- list.files(path = p, full.names = T)
  
  print(fn)
  
  df <- rbindlist(lapply(fn, f))
  
  p <- paste0("~/mr_results_snps/full/", n[y], "_mr_full.csv")
  
  print(p)
  
  write.csv(df, p, row.names = F, quote = F)
  
})
