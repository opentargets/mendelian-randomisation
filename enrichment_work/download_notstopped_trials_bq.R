# Using BigQuery (genomics-public-data:linkage_disequilibrium_1000G_phase_3.super_pop_EUR & open-targets-genetics:210608.variants)  for LD matrix construction

library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID

qry <- "SELECT DISTINCT 
    *
FROM
    `open-targets-genetics-dev.mohd_ds.davids_notstopped_annotated_ensid_efo`"

drugs <- bq_project_query(project, qry) %>% bq_table_download()
