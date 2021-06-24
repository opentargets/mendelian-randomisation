# extract table names from dataset to filter on

tbl <- bq_dataset_tables(paste(project_id, dataset_name2, sep = "."))

tblnames <- sapply(seq(1:length(tbl)), function (x) tbl[[x]]$table)

# # create bq dataset where filtered variants will exist
# 
# dataset_name3 <- paste0(dataset_name2, "_filtered")
# 
# if(bq_dataset_exists(paste0(project_id, ".", dataset_name3)) == FALSE) {
#   
#   ds <- bq_dataset(project_id, dataset_name3)
#   
#   bq_dataset_create(ds, location = "EU")
#   
# }

# dest_tbl <- paste(project_id, dataset_name3, tblnames, sep = ".")

filter_p <- function () {

sapply(seq(tblnames), function (x) {
  
  tbl <- tblnames[x]
  # dest_tbl <- dest_tbl[x]
  
  print(tbl)
  
  qry <- paste0("SELECT
  *
FROM
  `", dataset_name2, ".", tbl, "`
WHERE
  CAST(p as float64) <= 1e-5")
  
  # bq_perform_query(qry, 
  #                  billing = billing_project_id, 
  #                  destination_table = dest_tbl,
  #                  write_disposition = "WRITE_TRUNCATE",
  #                  priority = "BATCH")
  
  qry <- gsub("'", '"', qry)
  
  cmd <- paste0('bq query --replace --batch --destination_table ', dataset_name3, '.', tbl, ' --use_legacy_sql=false ',  "'", qry, "'")
  
  system(cmd)
  
})
  
}

