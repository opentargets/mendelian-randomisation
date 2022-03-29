library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID

qry1 <- "SELECT 
    study_id,
    anc
FROM 
    `open-targets-genetics-dev.mohd_ds.studies_efos` df1
FULL JOIN 
    `open-targets-genetics-dev.mohd_ds.parent_efos2` df2
ON 
    df1.trait_efos=df2.efo"

pefo <- bq_project_query(project, qry1) %>% bq_table_download(page_size = 30000) # toggle page_size if table is too large to be parsed

