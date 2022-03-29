library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID

qry1 <- "SELECT
    study_id,
    trait_reported,
    trait_category,
    trait_efos,
    n_initial,
    n_cases
FROM 
    `open-targets-genetics-dev.genetics_dev.studies`"

efo <- bq_project_query(project, qry1) %>% bq_table_download(page_size = 30000) # toggle page_size if table is too large to be parsed
