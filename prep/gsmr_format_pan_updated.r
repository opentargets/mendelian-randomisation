# LOAD LIBRARIES AND AUTHENTICATE

source("prep/loadLibraries_and_authenticate.r")

# LOAD CONFIG FILE

source("~/mr/mendelian-randomisation/configs/config_nightingale.r")


# CREATE DATASET IN BigQuery (BQ) PROJECT

dataset_name2 <- paste0(dataset_name, "_gsmr")

dataset_name3 <- paste0(dataset_name, "_gsmr", "_filtered")

if(bq_dataset_exists(paste0(project_id, ".", dataset_name)) == FALSE) {
  
  ds <- bq_dataset(project_id, dataset_name)
  
  bq_dataset_create(ds, location = "EU") 
  
  ds2 <- bq_dataset(project_id, dataset_name2)
  
  bq_dataset_create(ds2, location = "EU")
  
  ds3 <- bq_dataset(project_id, dataset_name3)
  
  bq_dataset_create(ds3, location = "EU")
  
}

cat("datasets created in BQ")


# LOAD TABLES IN BQ DATASET

# to check status of jobs > https://www.googleapis.com/bigquery/v2/projects/open-targets-ukbb/jobs/<job-id>

input1 <- strsplit(input_gc_bucket, "/")[[1]][3] # main bucket name

input2 <- gsub(paste("gs\\:\\/\\/", paste0(input1, "/"), sep = "|"), "", input_gc_bucket)

gc <- gcs_list_objects(bucket = input1, prefix = input2) # fetch file list from the buckets

cond <- nrow(bqr_list_tables(projectId = project_id, datasetId = dataset_name)) == nrow(gcs_list_objects(bucket = input1, prefix = input2))

if(cond==FALSE) {

if(parq==TRUE) {
  
  source("prep/load_pq_tables_bq.r")
  load_pq_tables_bq() # script loads tables from gc to the bq dataset
  
} else {
  
  source("prep/load_tables_bq.r")
  load_tables_bq() # script loads tables from gc to the bq dataset
  
}
  
}

cond <- nrow(bqr_list_tables(projectId = project_id, datasetId = dataset_name)) == nrow(gcs_list_objects(bucket = input1, prefix = input2))

if(cond==FALSE) {
  print("not all tables were loaded")
  break
}

cat("tables loaded in BQ")

# LOAD SCHEMA

if(schema=="autodetect") {
  
  p <- paste0(project_id, ":", dataset_name, ".", tbl2[1])
  
  system(paste0("bq show --format=prettyjson ", p, " | jq '.schema.fields' > ", "schema_", dataset_name))
  
  Sys.sleep(60) # since the above takes a bit of time after execcution
  
  schema <- jsonlite::read_json(paste0("schema_", dataset_name))
  
} else {
  
  schema <- try(readRDS(schema)) # schema (rds object)
  
  if(!exists("schema")) {
    
    if(grepl("parquet", schema_file)==TRUE) {
      
      schema_outcomes <- arrow::read_parquet(schema_file)
      
      schema <- schema_fields(schema_outcomes)
      
    } else {
      
      schema_outcomes <- data.table::fread(schema_file, stringsAsFactors = F)
      
      schema <- schema_fields(schema_outcomes)
      
    }
    
  }
  
}

# CREATE QUERY to (not necessarily in this order): 

# 1) deduplicate data based on variant ID, 
# 2) join on variant ID and annotate with UKB AF data (if use_ukb_eaf=TRUE), 
# 3) remove variants with null SE, 
# 4) re-annotate p-value = 0 with smallest double-precision float below which most statistical platforms assign null (10^-323 in BQ, python), # 5) select the right columns as appropriate

# bqt <- paste(project_id, dataset_name2, tbl2, sep = ".") # building BQ table references

cond <- nrow(bqr_list_tables(projectId = project_id, datasetId = dataset_name2)) == nrow(gcs_list_objects(bucket = input1, prefix = input2))

if(cond==FALSE) {

if(nchar(schema_effect_allele_freq)==0) {
  
  source("prep/process_tables_bq.r") # annotates eaf with ukb eaf
  process_tables_bq() # script process data in BQ to format and store in a separate bq dataset (<dataset>_gsmr)
  
} else {
  
  source("prep/process_tables_bq_noukbaf.r")
  process_tables_bq_noukbaf() # script process data in BQ to format and store in a separate bq dataset (<dataset>_gsmr)
  
}
  
}

cat("tables processed in BQ")

# P-VALUE FILTERED EXPOSURE TABLES

# source("loadLibraries_and_authenticate.r") # sometimes needs a re-authentication

cond <- nrow(bqr_list_tables(projectId = project_id, datasetId = dataset_name3)) == nrow(bqr_list_tables(projectId = project_id, datasetId = dataset_name))

if(cond==FALSE) {

source("prep/filter_p.r")

filter_p()

# Create a list of unique vars across all exposure GWASes

source("prep/create_unique_var_list.r")

create_unique_var_list()

# MATCHED OUTCOME TABLE

source("prep/curate_outcome.r")

curate_outcome()

}

# EXPORT

## EXPOSURES

source("prep/export_output_to_gcs_exposures.r")

export_output_to_gcs_exposures()

if(gb==TRUE) {
  
  gc_bucket <- output_gc_bucket
  
  source("prep/compose_and_delete_splits.r")
  
  compose_and_delete_splits()
  
}

## OUTCOMES

source("prep/export_output_to_gcs_outcomes.r")

export_output_to_gcs_outcomes()

if(gb==TRUE) {
  
  gc_bucket <- outcome_gc_bucket
  
  source("prep/compose_and_delete_splits.r")
  
  compose_and_delete_splits()
  
}
