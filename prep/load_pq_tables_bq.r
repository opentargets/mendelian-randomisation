# LOAD TABLES IN BQ DATASET

# FETCH FILE LIST FROM GOOGLE CLOUD (GC) BUCKETS

input1 <- strsplit(input_gc_bucket, "/")[[1]][3] # main bucket name

input2 <- gsub(paste("gs\\:\\/\\/", paste0(input1, "/"), sep = "|"), "", input_gc_bucket)

gc <- gcs_list_objects(bucket = input1, prefix = input2) # fetch file list from the buckets

if(grepl(".parquet", gc$name)[1]==TRUE) {
  
  gc$url <- sapply(seq(gc$name), function (x) paste0("gs://", input1, "/", input2, "/", strsplit(gc$name[x], "/")[[1]][2]))
  gc <- gc[!duplicated(gc$url),]
  
} else {
  
  gc$url <- paste0('gs://', input1, '/', gc$name) # annotate with GC URLs
  
}

# CREATE TABLE NAMES (in BQ and for GCS OUTPUT)

l <- length(strsplit(gc$url[1], "/")[[1]])

l2 <- length(strsplit(strsplit(gc$url[1], "/")[[1]][l], "\\.")[[1]])

if (format_suffix==TRUE) {
  
  format_suffix <- paste0(".", strsplit(strsplit(gc$url[1], "/")[[1]][l], "\\.")[[1]][[l2]])
  
  tbl <- sapply(seq(gc$url), function (x) gsub(format_suffix, "", strsplit(gc$url[x], "/")[[1]][l])) # generating tbl names
  
} else {
  
  tbl <- sapply(seq(gc$url), function (x) strsplit(gc$url[x], "/")[[1]][l])  # generating tbl names
  
}

tbl2 <- gsub("\\.", "_", tbl) # if filenames have ".", will change this to "underscore", other name will collide with the way BQ references project.dataset.table names

bqt <- paste(project_id, dataset_name, tbl2, sep = ".")

load_pq_tables_bq <- function () {
  
  require(bigrquery)
  # bq_auth(path = json_path)
  
  sapply(seq(gc$url), function (x) { # loop to load tables in BQ project dataset from GC
    
    url <- gc$url[x]
    url <- paste0(gc$url[x], "/*.parquet")
    # tbl <- tbl2[x]
    
    print(url)
    #print(tbl)
    
    bq_perform_load(bqt[x],
                    url,
                    billing = billing,
                    source_format = "PARQUET",
                    write_disposition = "WRITE_TRUNCATE")
    
    # bqr_upload_data(projectId = project,
    #                 datasetId = ds_name,
    #                 tableId = tbl,
    #                 upload_data = url,
    #                 writeDisposition = "WRITE_TRUNCATE",
    #                 create = "CREATE_IF_NEEDED",
    #                 sourceFormat = "CSV")
    
    # print("table loaded from GCS to BQ")
    
  })
  
}

# tables <- bqr_list_tables(projectId = project, datasetId = ds_name) # verify if tables have loaded (the table loading might take some time)

# sapply(seq(gc$url)[1:2], function (x) { # loop to load tables in BQ project dataset from GC
#   
#   url <- gc$url[x]
#   
#   print(gc$url[x])
#   
#   bq_perform_load(bqt[x], url, billing = billing, source_format = "CSV", write_disposition = "WRITE_TRUNCATE", fields = schema, nskip = 1) # auto-detect schema. if schema specified, for some reason, it cannot load table and generates errors
#   
#   # print("table loaded from GCS to BQ")
#   
# }) # this uses bigRquery, but omitted because it does not have argument to specify delimiter