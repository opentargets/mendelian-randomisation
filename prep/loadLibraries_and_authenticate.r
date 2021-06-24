# Load libraries

library(bigrquery)
library(bigQueryR)
library(dplyr)
library(data.table)
library(googleCloudStorageR)
httr::set_config(httr::config(http_version = 0))

# authenticate

json_path = "open-targets-genetics-fc5b6cda58e5.json"

bq_auth(path = json_path, use_oob = T)
bqr_auth(json_file = json_path)
gcs_auth(json_file = json_path)