# Creating table names

dataset_name3 <- paste0(dataset_name2, "_filtered")

dest <- paste0(dataset_name3)

qry1 <- paste0("SELECT * FROM `", dest, ".__TABLES__`")

tbl2 <- bq_project_query(project_id, qry1) %>% bq_table_download() %>% select(table_id, size_bytes)

gb <- any(ifelse(tbl2$size_bytes > 1e9, "TRUE", "FALSE")==TRUE) # this will determine whether a script with wildcards is used, as BQ splits while exporting to GC any file  > 1 GB

# export function

export_output_to_gcs_exposures <- function () {
  
  if(gb==FALSE) {
    
    sapply(seq(tbl2$table_id), function (x) {
      
      # require(bigrquery)
      
      tbl <- tbl2$table_id[x]
      
      print(tbl)
      
      # tbl_size <- bq_table_size(paste(project, ds_name2, tbl, sep = ".")) # wanted to specify if tables size is below 1GB, export as single file rather than split with trailing zeros, but the export option can only exporrt with wildcards. So will have to compose and delete splits later on.
      
      #   qry2 <- paste0("EXPORT DATA
      #   OPTIONS(uri='", output_gc_bucket, tbl, "*.csv', #files bigger than 1 GB will be split by BQ as it sends the processed files to GC
      #     format='CSV',
      #     overwrite=TRUE,
      #     header=FALSE,
      #     field_delimiter=' ') AS
      # SELECT
      #   *
      # FROM
      #   `", dataset_name2, ".", tbl, "`")
      #   
      #   bq_perform_query(qry2, 
      #                    billing = billing_project, 
      #                    priority = "BATCH")
      
      
      origin_tbls <- paste0(project_id, ":", dest, ".", tbl)
      
      dest_tbls <- paste0(output_gc_bucket, tbl, ".csv")
      
      cmd <- paste0('bq --location=EU extract --destination_format CSV --compression NONE --field_delimiter " " --print_header=false ', origin_tbls,' ', dest_tbls)
      
      system(cmd)
      
      
      # print("updated tables sent to GCS")
      
    })
    
  } else {
    
    sapply(seq(tbl2$table_id), function (x) {
      
    #  require(bigrquery)
      
      tbl <- tbl2$table_id[x]
      
      print(tbl)
      
      origin_tbls <- paste0(project_id, ":", dest, ".", tbl)
      
      dest_tbls <- paste0(output_gc_bucket, tbl, "*.csv") # using wildcards to export > 1GB files
      
      cmd <- paste0('bq --location=EU extract --destination_format CSV --compression NONE --field_delimiter " " --print_header=false ', origin_tbls,' ', dest_tbls)
      
      system(cmd)
      
      # print("updated tables sent to GCS")
      
    })
    
    
  }
  
}