


compose_and_delete_splits <- function ()  {

output1 <- strsplit(gc_bucket, "/")[[1]][3] # main bucket name

output2 <- gsub(paste("gs\\:\\/\\/", paste0(output1, "/"), sep = "|"), "", gc_bucket)

# buckets <- gcs_list_buckets(project) # fetch all buckets from projects

gc_output <- gcs_list_objects(bucket = output1, prefix = output2)

l <- length(strsplit(gc_output$name[1], "/")[[1]])

gc_output$table <- sapply(seq(gc_output$name), function (x) gsub(".csv", "", strsplit(gc_output$name[x], "/")[[1]][l]))

gc2 <- gc_output[grep("00000000000", gc_output$name),]

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
} # extract last 12 characters

gc2$dfn <- sapply(seq(gc2$table), function (x) gsub(substrRight(gc2$table[x], 12), "", gc2$table[x]))

length(unique(gc2$dfn))==length(tbl2)-1 # for outcomes,"GCST004692" study was not present in gs://genetics-portal-sumstats-b38/unfiltered/gwas/, so did not copy in gs://genetics-portal-ukbb-mr-eur/outcomes/, so expect to have length(tbl2)-1

implode <- function(..., sep=' ') {
  paste(..., collapse=sep)
}

gc2$dfn_comp <- sapply(seq(gc2$dfn), function (x) {

  obj <-gc2$dfn[x]

  print(obj)

  eg <- gc2[gc2$dfn==obj,]
  eg2 <- paste(eg$table, sep = " ")
  vec1 <- paste0(gc_bucket, eg2, ".csv")

  vec2 <- implode(vec1)
  vec3 <- paste0(vec2, " ", gc_bucket, obj, ".csv")


})

gc3 <- unique(gc2$dfn_comp) # vector of datasets to compose with their correct names attached

write.table(gc3, "list_to_recompose.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

gc4 <- paste0(gc_bucket, gc2$table, ".csv")

write.table(gc4, "split_files_to_delete.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

system("chmod +x ./compose_only.sh && ./compose_only.sh")

system("chmod +x ./delete_splits.sh && ./delete_splits.sh")

}