# Using BigQuery (genomics-public-data:linkage_disequilibrium_1000G_phase_3.super_pop_EUR & open-targets-genetics:210608.variants)  for LD matrix construction

library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID

# qry <- "SELECT DISTINCT
#     *
# FROM
#     `open-targets-genetics-dev.mohd_ds.davids_stopped_annotated_ensid_efo`"

qry <- "# Aim: To link NCT IDs from David's stopped dataset with ensid-efo trait pairs linked via chembl IDs

WITH bqdf2 AS (
WITH bqdf AS (
SELECT DISTINCT # annotating with chembl_id and disease efo
    df1.*,
    df2.chembl_id,
    df2.diseaseFromSource,
    df2.diseaseFromSourceMappedId
FROM
    `open-targets-genetics-dev.mohd_ds.stopped_predictions06_02` df1
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
    bqdf2.*,
    df4.parent_efo
FROM
    bqdf2
LEFT JOIN
    `open-targets-genetics-dev.mohd_ds.parent_efos2` df4
ON
    bqdf2.diseaseFromSourceMappedId=df4.direct_efo"

drugs <- bq_project_query(project, qry) %>% bq_table_download()
