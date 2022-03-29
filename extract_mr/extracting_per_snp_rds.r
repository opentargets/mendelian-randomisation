# Load library and authenticate

# source("loadLibraries_and_authenticate.r")

library(googleCloudStorageR)

gcs_auth(json_file = "open-targets-genetics-fc5b6cda58e5.json")

# Process the gz file (this takes a while!) but needs to be done only once (the resulting rds file loads faster)

# WHERE PROTEIN/METABOLITES ARE NOT IN MULTIPLE DIRECTORIES (i.e. except Sun and Nightingale)

prot <- c("KETTUNEN", "PIETZNER", "SCALLOP", "SUHRE", "FOLKERSEN", "HILLARY", "OLLI", "DRAISMA", "SHIN")

# prot <- "KETTUNEN"

# prot <- "SUN"
# 
# a <- seq(1,11,1)
# 
# dir <- sprintf("dir_%03d", a)

direction <- "forward"

ver <- "v2"

# dirs <- gcs_list_objects(bucket = "genetics-portal-ukbb-mr-eur", prefix = paste0(direction, "-pan-gsmr-outputs/", prot))

source("gsmr_plot.r") # hacked the gsmr_snp_effect function to show MR results with less than 3 instruments

parallel::mclapply(seq(prot), function (x) {

  system(paste0("mkdir ", prot[x], "_gsmr_",ver, "_data"))

  gcs_get_object(paste0("gs://genetics-portal-ukbb-mr-eur/", direction, "-pan-gsmr-outputs-nsnp1_", ver, "/", prot[x], "/results.eff_plot.gz"), saveToDisk = paste0(prot[x], "_gsmr_",ver, "_data/", prot[x], "_results_gsmr_snp.gz"))

  mrdf <- read_gsmr_data(paste0(prot[x], "_gsmr_", ver, "_data/", prot[x], "_results_gsmr_snp.gz")) # this takes a bit of time

  saveRDS(mrdf, paste0(prot[x], "_gsmr_",ver, "_data/", prot[x], "_results_gsmr_snp.rds"))

}, mc.cores = parallel::detectCores())

# SUN and NIGHTINGALE (where their proteins are in multiple directories)

prot2 <- c("SUN", "NIGHTINGALE")

ver <- "v2"

dir <- sprintf("dir_%03d", seq(1,11,1)) # 11 directories

lapply(seq(prot2), function (y) {
  
  require(googleCloudStorageR)
  
  gcs_auth(json_file = "open-targets-genetics-fc5b6cda58e5.json")
  
  prot <- prot2[y]
  
  print(prot)
  
  # system(paste0("mkdir ", prot, "_gsmr_v2_data"))
  
  parallel::mclapply(seq(dir), function (x) {
    
    gcs_get_object(paste0("gs://genetics-portal-ukbb-mr-eur/", direction, "-pan-gsmr-outputs-nsnp1_", ver, "/", prot, "/", dir[x], "/results.eff_plot.gz"), saveToDisk = paste0(prot, "_gsmr_",ver, "_data/", dir[x], "_results_gsmr_snp.gz"))
    
    mrdf <- read_gsmr_data(paste0(prot, "_gsmr_",ver, "_data/", dir[x], "_results_gsmr_snp.gz")) # this takes a bit of time
    
    saveRDS(mrdf, paste0(prot, "_gsmr_",ver, "_data/", dir[x], "_results_gsmr_snp.rds"))
    
  }, mc.cores = parallel::detectCores())
  
})

# special run for NIGHTINGALE v1 

# prot <- "N"

# system(paste0("mkdir ", prot[4], "_gsmr_data"))
# #
# gcs_get_object(paste0("gs://genetics-portal-ukbb-mr-eur/", direction, "-pan-gsmr-outputs-nsnp1/", prot[4], "/results.eff_plot.gz"), saveToDisk = paste0(prot[4], "_gsmr_data/", prot[4], "_results_gsmr_snp.gz"))
# 
# mrdf <- read_gsmr_data(paste0(prot[4], "_gsmr_data/", prot[4], "_results_gsmr_snp.gz")) # this takes a bit of time
# 
# saveRDS(mrdf, paste0("~/", prot[4], "_gsmr_data/", prot[4], "_results_gsmr_snp.rds"))

