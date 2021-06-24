out_tblnames <- bqr_list_tables(projectId = project_id, datasetId = outcome_ds)

out_tblnames <- out_tblnames[-grep("parquet", out_tblnames$tableId),] # remove files labelled as parquet (this is a non-generalizable step applied only for outcomes_gsmr dataset)

# ds2 <- bq_dataset(project_id, new_outcome_ds)
# 
# bq_dataset_create(ds2, location = "EU")
# 
# out_tbls <- bq_dataset_tables(paste(project_id, outcome_ds, sep = "."))

# out_tblnames <- sapply(seq(1:length(out_tbls)), function (x) out_tbls[[x]]$table)

out_tblnames <- out_tblnames$tableId

unique_var_tbl <- paste0(unique_var_ds, ".", dataset_name3, "_unique_snps")

dest <- paste0("outcomes", "_", dataset_name3)

if(bq_dataset_exists(paste0(project_id, ".", dest)) == FALSE) {
  
  ds <- bq_dataset(project_id, dest)
  
  bq_dataset_create(ds, location = "EU")
  
}

# out_dest_tbl <- paste(project_id, new_outcome_ds, out_tblnames, sep = ".")

curate_outcome <- function () {

sapply(seq(out_tblnames), function (x) {
  
  tbl <- out_tblnames[x]
  # dest_tbl <- out_dest_tbl[x]
  
  print(tbl)
  
  qry <- paste0("WITH
  ds AS (
  SELECT
    *
  FROM
    `", outcome_ds, ".", tbl, "` m
  INNER JOIN
    `", unique_var_tbl, "` g
  ON
    m.SNP = g.VAR)
SELECT
  * EXCEPT (VAR)
FROM
  ds")
  
  # bq_perform_query(qry, 
  #                  billing = billing_project_id, 
  #                  destination_table = dest_tbl,
  #                  write_disposition = "WRITE_TRUNCATE",
  #                  priority = "BATCH")
  
  qry <- gsub("'", '"', qry)
  
  cmd <- paste0('bq query --replace --batch --destination_table ', dest, '.', tbl, ' --use_legacy_sql=false ',  "'", qry, "'")
  
  system(cmd)
  
})
  
}
