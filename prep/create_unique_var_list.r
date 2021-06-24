
origin_tbls <- paste(project_id, dataset_name3, "*", sep = ".")

tbl <- paste0(dataset_name3, "_unique_snps")

# dest_tbl <- paste(project_id, unique_var_ds, paste0(dataset_name3, "_unique_snps"), sep = ".")

create_unique_var_list <- function () {

qry <- paste0("WITH
  TBL AS (
SELECT
  SNP as VAR,
  ROW_NUMBER() OVER (PARTITION BY SNP) row_number
FROM
  `", origin_tbls, "`)
SELECT
  VAR
FROM
  TBL
WHERE
  row_number = 1")

qry <- gsub("'", '"', qry)

cmd <- paste0('bq query --batch --destination_table ', unique_var_ds, '.', tbl, ' --use_legacy_sql=false ',  "'", qry, "'")

system(cmd)

}
