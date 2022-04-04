# Using BigQuery (genomics-public-data:linkage_disequilibrium_1000G_phase_3.super_pop_EUR & open-targets-genetics:210608.variants)  for LD matrix construction

library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID

qry <- paste0(c("WITH bqdf2 AS (
WITH bqdf AS (
SELECT DISTINCT # annotating with chembl_id and disease efo
    df1.*,
    df2.chembl_id,
    df2.diseaseFromSourceMappedId
FROM 
    `open-targets-genetics-dev.mohd_ds.notstopped_predictions06_02` df1
LEFT JOIN 
    `open-targets-genetics-dev.mohd_ds.chembl_NCTs` df2
ON
    df1.CT=df2.CT
WHERE 
    df2.chembl_id IS NOT NULL
)
SELECT DISTINCT # annotating with ensid
    bqdf.*,
    df3.ensid
FROM 
    bqdf 
LEFT JOIN 
    `open-targets-genetics-dev.mohd_ds.chembl_id_ensid` df3
ON
    bqdf.chembl_id=df3.chembl_id
)
SELECT DISTINCT # annotating with parent efos
    bqdf2.CT,
    bqdf2.diseaseFromSourceMappedId,
    df4.trait_reported,
    df4.trait_category
FROM 
    bqdf2 
LEFT JOIN 
    `open-targets-genetics-dev.mohd_ds.parent_efos2` df4
ON
    bqdf2.diseaseFromSourceMappedId=df4.direct_efo
WHERE 
    bqdf2.CT IN (", toString(shQuote(ncts)),")"), collapse = "") 

nct_traits <- bq_project_query(project, qry) %>% bq_table_download()
