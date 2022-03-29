

# linking trait name to hmdb identifiers

# grep -oP "<name>Alanine</name>" serum_metabolites.xml | cut -d ">" -f 2 | cut -d "<" -f 1

library("xml2")
library(purrr)
# library(furrr)
# future::plan(multiprocess)

x <- read_xml("~/hmdb/serum_metabolites.xml")

# x_list <- as_list(x)

# xdf <- x %>% xml_find_all('.//d1:metabolite') %>% future_map(as_list) %>% future_map_dfc(flatten)

# xdf <- x %>% xml_find_all('.//d1:metabolite') %>% map(as_list) %>% map_df(flatten)


#get metabolite nodes (only first three used in this sample)

met.nodes <- xml_find_all(x, ".//d1:metabolite")

# metlst <- as_list(met.nodes)

#list of data.frames with secondary accessions, synonyms, and uniprot ID

xpath_child.v <- c("./d1:secondary_accessions/d1:accession",
                   "./d1:name",
                    "./d1:synonyms/d1:synonym",
                   "./d1:protein_associations/d1:protein/d1:uniprot_id")

#what names should they get in the list?
child.names.v <- c("secondary_accessions",
                   "name",
                    "synonyms",
                   "uniprot")

# temp <- lapply( xpath_child.v, function(y) { 
#   xml_find_all(met.nodes[10000], y ) %>% xml_text() %>% data.frame(value = ., stringsAsFactors = F)
# })
# 
# names(temp) = child.names.v
# 
# eg <- xml_find_all(met.nodes[10000], "./d1:protein_associations/d1:protein/d1:uniprot_id") %>% xml_text() %>% data.frame(value = ., stringsAsFactors = F)

#first, loop over the met.nodes (inspired from https://stackoverflow.com/questions/63392837/convert-xml-file-to-a-list-of-data-frames-in-r)

dflst <- parallel::mclapply(met.nodes, function(x) {
  
  #second, loop over the xpath desired child-nodes
  
  temp <- lapply( xpath_child.v, function(y) { 
    xml_find_all(x, y) %>% xml_text()
  })
  
  #set their names
  
  # temp <- bind_rows(temp)
  
  # temp <- data.frame(temp, stringsAsFactors = F)
  
  names(temp) = child.names.v
  
  # temp <- as.data.frame(temp, stringsAsFactors = FALSE)
  
  return(temp)
  
}, mc.cores = parallel::detectCores())

saveRDS(dflst, "~/hmdb_serum.rds")

dflst2 <- parallel::mclapply(seq(dflst), function (x) {
  
  df <- dflst[[x]]
  
  df2 <- data.frame(t(plyr::ldply(df, rbind)), stringsAsFactors = F)
  
  names(df2) <- df2[1,]
  
  df3 <- df2[-1,]
  
  df3$name <- df3$name[1]
  
  df3
  
}, mc.cores = parallel::detectCores())


df <- bind_rows(dflst2)


saveRDS(df, "hmdb_serum_df.rds")


