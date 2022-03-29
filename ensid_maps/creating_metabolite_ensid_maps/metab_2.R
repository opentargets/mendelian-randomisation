# Using the rds from metab.R to annotate metabolites from metabolite GWAS

library(dplyr)

# Loading "hmdb_serum_df.rds"

hmdb <- readRDS("~/hmdb_serum_df.rds")

# Generating metabolite data in the same format as protein data

# metabs <- c("SHIN_2014", "KETTUNEN_2016", "DRAISMA_2015")

# metab_files <- grep(paste(metabs, collapse = "|"), dir("~/mr_results_snps/2_annotation_eggerwm_outcometraits/", full.names = T), value = T)

# dfm <- readRDS(metab_files[3])

# SHIN map

ws <- rio::import("https://static-content.springer.com/esm/art%3A10.1038%2Fng.2982/MediaObjects/41588_2014_BFng2982_MOESM50_ESM.xlsx", readxl = F, which = 3) %>% select(1:2)

ws <- ws[-c(1:2),]

names(ws) <- c("exposure", "trait")

# remove asterisk and parenthesis which can interfere with pattern matching

ws$trait <- gsub("\\*", "", ws$trait)

ws$trait <- gsub("\\s*\\([^\\)]+\\)\\s*$", "", ws$trait) # only remove parenthesis and words within that are the end of the trait strings

df <- parallel::mclapply(seq(ws$trait), function (x) {
  
  pos <- try(grep(paste0("^", ws$trait[x], "$"), hmdb$synonyms, ignore.case = TRUE))
  
  eg1 <- try(hmdb[hmdb$name==hmdb$name[pos],])
  
  # eg1 <- try(hmdb[grep(paste0("^", ws$trait[1], "$"), hmdb$synonyms, ignore.case = TRUE),])
  
  eg2 <- try(hmdb[grep(paste0("^", ws$trait[x], "$"), hmdb$name, ignore.case = TRUE),])
  
  eg3 <- rbind(eg1, eg2)
  
  eg <- eg3
  
  # eg <- eg3[!duplicated(eg3),]
  
  ws2 <- ws[ws$trait==ws$trait[x],]
  
  if(is(eg, "try-error")) {
    
    eg <- data.frame(secondary_accessions = NA, name = NA, synonyms = NA, uniprot = NA, exposure = ws2$exposure[1], trait = ws2$trait[1], Data = "SHIN_2014", stringsAsFactors = F)
  
  } else if(nrow(eg)!=0) {
    
    eg <- merge(eg, ws2, all = TRUE)
    eg$Data <- "SHIN_2014"
    
    eg <- data.frame(eg, stringsAsFactors = F)
  
  } else {
    
  eg <- data.frame(secondary_accessions = NA, name = NA, synonyms = NA, uniprot = NA, exposure = ws2$exposure[1], trait = ws2$trait[1], Data = "SHIN_2014", stringsAsFactors = F)
    
  }
  
  
}, mc.cores = parallel::detectCores(), mc.preschedule = FALSE)

# df2 <- do.call(rbind, df)

df2 <- bind_rows(df)

#df2$V1 <- NULL

df2$source <- with(df2, ifelse(is.na(uniprot), NA, "hmdb"))

# df4 <- df3 %>% select(Data, exposure, trait, secondary_accessions, name, synonyms, uniprot)
# 
# names(df4) <- c("metab_data", "metab_exposure_id", "metab_exposure_trait", "hmdb_accession", "hmdb_name", "hmdb_synonyms", "hmdb_uniprot")

# row.names(df2) <- 1:nrow(df2)

saveRDS(df2, "shin_map_hmdb.rds")

# manual mapping 

# mets <- unique(df3$exposure) # 224 out of 529 metabolites

# unmapped_df <- ws[!ws$exposure %in% mets,] # manually mapping these or post-MR mapping (if they are 5% FDR significant)

# The large majority of the unmapped in Shin are X-* metabolites + a few where hmdb was not able to assign a gene. Will not be taking this for GO annotations as there is high chance of unreliable annotations, but will deal with annotations post-MR analysis on an ad-hoc basis.

# pos_overlap <- lapply(seq(unmapped_df$trait3), function (x) try(grep(unmapped_df$trait3[x], df3$name, ignore.case = T)))
# 
# pos_overlap2 <- which(length(pos_overlap)!=0)

# N-acetylalanine

# creatine

# creatinine

# pyroglutamine

# cysteine-glutathione disulfide


# KETTUNEN map

# url pdf (page 161 - 162) > https://static-content.springer.com/esm/art%3A10.1038%2Fncomms11122/MediaObjects/41467_2016_BFncomms11122_MOESM71_ESM.pdf (extracted pages 161-162 and converted to excel online)

ket <- rio::import("~/kettunen_map.xlsx", readxl = FALSE, colNames = FALSE, startRow = 2)

names(ket) <- ket[1,]

ket <- ket[-c(1, 125),]

ws <- ket[,c(1,2)]

names(ws) <- c("exposure", "trait")

# remove asterisk and parenthesis which can interfere with pattern matching

ws$trait <- gsub("\\*|\\s*\\([^\\)]+\\)\\s*$", "", ws$trait) # only remove parenthesis and words within that are the end of the trait strings

# df <- parallel::mclapply(seq(ws$trait2), function (x) {
#   
#   eg1 <- try(hmdb[grep(paste0("^", ws$trait2[x], "$"), hmdb$synonyms, ignore.case = TRUE),])
#   
#   eg2 <- try(hmdb[grep(paste0("^", ws$trait2[x], "$"), hmdb$name, ignore.case = TRUE),])
#   
#   eg3 <- rbind(eg1, eg2)
#   
#   eg <- eg3[!duplicated(eg3),]
#   
#   ws2 <- ws[ws$trait2==ws$trait2[x],]
#   
#   if(nrow(eg)!=0) {
#     
#     eg$exposure <- ws2$exposure[1]
#     eg$trait <- ws2$trait3[1]
#     eg$Data <- "KETTUNEN_2016"
#     
#     eg <- data.frame(eg, stringsAsFactors = F)
#     
#   } else {
#     
#     eg <- data.frame(secondary_accessions = NA, name = NA, synonyms = NA, uniprot = NA, exposure = ws2$exposure[1], trait = ws2$trait3[1], Data = "KETTUNEN_2016", stringsAsFactors = F)
#     
#   }
#   
#   eg
#   
# }, mc.cores = parallel::detectCores(), mc.preschedule = FALSE)
# 
# df2 <- do.call(rbind, df)
# 
# df2$source <- with(df2, ifelse(is.na(df2$uniprot), NA, "hmdb"))

df <- parallel::mclapply(seq(ws$trait), function (x) {
  
  pos <- try(grep(paste0("^", ws$trait[x], "$"), hmdb$synonyms, ignore.case = TRUE))
  
  eg1 <- try(hmdb[hmdb$name==hmdb$name[pos],])
  
  eg2 <- try(hmdb[grep(paste0("^", ws$trait[x], "$"), hmdb$name, ignore.case = TRUE),])
  
  eg3 <- rbind(eg1, eg2)
  
  eg <- eg3
  
  # eg <- eg3[!duplicated(eg3),]
  
  ws2 <- ws[ws$trait==ws$trait[x],]
  
  if(is(eg, "try-error")) {
    
    eg <- data.frame(secondary_accessions = NA, name = NA, synonyms = NA, uniprot = NA, exposure = ws2$exposure[1], trait = ws2$trait[1], Data = "KETTUNEN_2016", stringsAsFactors = F)
    
  } else if(nrow(eg)!=0) {
    
    eg <- merge(eg, ws2, all = TRUE)
    eg$Data <- "KETTUNEN_2016"
    
    eg <- data.frame(eg, stringsAsFactors = F)
    
  } else {
    
    eg <- data.frame(secondary_accessions = NA, name = NA, synonyms = NA, uniprot = NA, exposure = ws2$exposure[1], trait = ws2$trait[1], Data = "KETTUNEN_2016", stringsAsFactors = F)
    
  }
  
  
}, mc.cores = parallel::detectCores(), mc.preschedule = FALSE)

# df2 <- do.call(rbind, df)

df2 <- bind_rows(df)

# df2$V1 <- NULL

df2$source <- with(df2, ifelse(is.na(uniprot), NA, "hmdb"))


# df3 <- df2[-grep("Error in", df2$exposure),]

mets1 <- unique(df2$exposure[which(is.na(df2$source))])

mets2 <- unique(df2$exposure[which(df2$source=="hmdb")])

mets3 <- setdiff(mets1, mets2)

unmapped_df <- ws[ws$exposure %in% mets3,]

# manual mapping a few

row.names(unmapped_df) <- 1:nrow(unmapped_df)

unmapped_df$uniprot <- ""

unmapped_df$uniprot <- as.list(unmapped_df$uniprot) # will expand rows later

unmapped_df$uniprot[grep("^Alb$", unmapped_df$exposure)] <- "P02768"

unmapped_df$uniprot[grep("^ApoA1$", unmapped_df$exposure)] <- "P02647"

unmapped_df$uniprot[grep("^ApoB$", unmapped_df$exposure)] <- "P04114"

unmapped_df$uniprot[grep("ratio", unmapped_df$exposure)] <- "ratio" # labelling to remove this later on

unmapped_df$uniprot[grep("^DHA$", unmapped_df$exposure)] <- list(c("O00154", "P49753", "Q8N9L9", "Q7Z449", "Q86TX2")) # https://hmdb.ca/metabolites/HMDB0002183

unmapped_df$uniprot[grep("^bOHBut$", unmapped_df$exposure)] <- list(c("O00204", "P22309", "Q02338", "Q9BUT1", "Q9NT06")) # https://hmdb.ca/metabolites/HMDB0000011

# using QuickGo API and evidence codes: ECO:0000269, ECO:0006056, ECO:0000250, ECO:0000245, ECO:0000303, ECO:0000305 (http://wiki.geneontology.org/index.php/Guide_to_GO_Evidence_Codes)

mets <- c("^Tot.FA$", "^CH2.in.FA$", "^DB.in.FA$", "^FALen$", "^FAw3$", "^FAw6$", "^FAw79S$", "^TotPG$", "^Gp$", "HDL", "IDL", "LDL", "VLDL", "^LA$", "^Est.C$", "^Free.C$", "^Serum.C$", "^MUFA$", "^otPUFA$", "^PC$", "^Serum.TG$", "^SM$", "^Crea$", "^Lac$", "^Leu$", "^Val$") 

mets_go <- c("GO:0006631", "GO:0006631", "GO:0006631", "GO:0006631", "GO:0006631", "GO:0006631", "GO:0006631", "GO:0006506", "GO:0006473", "GO:0034364", "GO:0034363", "GO:0034362", "GO:0034361", "GO:0043651", "GO:0034435", "GO:0008203", "GO:0008203", "GO:1903964", "GO:0034626", "GO:0006656", "GO:0006641", "GO:0006684", "GO:0046449", "GO:0006089", "GO:0006551", "GO:0006573")

metsdf <- data.frame(metabolites = mets, goterms = mets_go, stringsAsFactors = F)

library(httr)
library(jsonlite)
library(xml2)

lapply(seq(metsdf$metabolites), function (x) {
  
  go <- gsub(":", "%3A", metsdf$goterms[x])
  
  print(go)
  
  requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/annotation/search?goId=", go, "&taxonId=9606&evidenceCode=ECO%3A0000269%2C%20ECO%3A0006056%2C%20ECO%3A0000250%2C%20ECO%3A0000245%2C%20ECO%3A0000303%2C%20ECO%3A0000305&limit=200")
  
  r <- GET(requestURL, accept("application/json"))
  
  # stop_for_status(r)
  
  json <- toJSON(content(r))
  
  eg <- fromJSON(json)
  
  eg2 <- eg$results
  
  eg3 <- as.character(unique(eg2$geneProductId))
  
  eg4 <- gsub('.*:(.*)','\\1', eg3)
  
  # metsdf$uniprot[x] <<- list(eg4)
  
  unmapped_df$uniprot[grep(metsdf$metabolites[x], unmapped_df$exposure)] <<- list(eg4)
  
})

df3 <- unmapped_df[-which(unmapped_df$uniprot=="ratio"),]

df3$secondary_accessions <- NA
df3$name <- NA
df3$synonyms <- NA
df3$Data <- "KETTUNEN_2016"
df3$source <- with(df3, ifelse(is.na(df3$uniprot), NA, "GO"))

df4 <- df3 %>% select(secondary_accessions, name, synonyms, uniprot, exposure, trait, Data, source)

# final df

df <- rbind(df2, df4)

# group by exposure trait, and ret

saveRDS(df, "kettunen_map_hmdb_go.rds")

# go <- "GO%3A0034364"
# 
# requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/annotation/search?goId=", go, "&taxonId=9606&evidenceCode=ECO%3A0000269%2C%20ECO%3A0006056%2C%20ECO%3A0000250%2C%20ECO%3A0000245%2C%20ECO%3A0000303%2C%20ECO%3A0000305&limit=200")
# 
# r <- GET(requestURL, accept("application/json"))
# 
# stop_for_status(r)
# 
# json <- toJSON(content(r))
# 
# eg <- fromJSON(json)
# 
# eg2 <- eg$results
# 
# eg3 <- as.character(unique(eg2$geneProductId))
# 
# eg4 <- gsub('.*:(.*)','\\1', eg3)
# 
# # fatty acid metabolism > https://www.ebi.ac.uk/QuickGO/annotations?goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in&goId=GO:0006631&taxonId=9606&taxonUsage=exact
# 
# fa <- data.table::fread("~/QuickGO-annotations-1631803045463-20210916.tsv", stringsAsFactors = F)
# 
# fagenes <- unique(fa$`GENE PRODUCT ID`) # 632 genes
# 
# unmapped_df$uniprot[grep("fatty acid", unmapped_df$trait2, ignore.case = T)] <- list(fagenes)
# 
# unmapped_df$uniprot[grep("^Crea$", unmapped_df$exposure)] <- list(c("Q92481", "X6RBG4")) # https://www.ebi.ac.uk/QuickGO/annotations?goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in&goId=GO:0097273&taxonId=9606&taxonUsage=descendants
# 
# 
# # protein acetylation > https://www.ebi.ac.uk/QuickGO/annotations?goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in&goId=GO:0006473&taxonId=9606&taxonUsage=descendants
# 
# glycoprotein_acetyl <- data.table::fread("~/QuickGO-annotations-1631806000445-20210916.tsv", stringsAsFactors = F)
# 
# gagenes <- unique(glycoprotein_acetyl$`GENE PRODUCT ID`) # 343 genes
# 
# unmapped_df$uniprot[grep("^Gp$", unmapped_df$exposure)] <- list(gagenes)
# 
# # HDL (doing the annotations more systematically from here) (GO:0034364)
# 
# library(httr)
# library(jsonlite)
# library(xml2)
# 
# go <- "GO%3A0034364"
# 
# requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/annotation/search?goId=", go, "&taxonId=9606&evidenceCode=ECO%3A0000269%2C%20ECO%3A0006056%2C%20ECO%3A0000250%2C%20ECO%3A0000245%2C%20ECO%3A0000303%2C%20ECO%3A0000305&limit=200")
# 
# r <- GET(requestURL, accept("application/json"))
# 
# stop_for_status(r)
# 
# json <- toJSON(content(r))
# 
# eg <- fromJSON(json)
# 
# eg2 <- eg$results
# 
# eg3 <- as.character(unique(eg2$geneProductId))
# 
# eg4 <- gsub('.*:(.*)','\\1', eg3)
# 
# unmapped_df$uniprot[grep("HDL", unmapped_df$exposure)] <- list(eg4)
# 
# # IDL
# 
# # LDL
# 
# # VLDL
# 
# # phosphoglycerides
# 
# # Urea
# 
# # sphingomyelins
# 
# # Phosphatidylcholine
# 
# # linoleic acid
# 
# # triglycerides
# 
# # total cholesterol & esterified cholesterol

  
# DRAISMA map

ws <- rio::import("~/draisma_map_file.xlsx") %>% select(1,2)

names(ws) <- c("exposure", "trait")

# remove asterisk and parenthesis which can interfere with pattern matching

ws$trait <- gsub("\\*|\\s*\\([^\\)]+\\)\\s*$", "", ws$trait) # only remove parenthesis and words within that are the end of the trait strings


df <- parallel::mclapply(seq(ws$trait), function (x) {
  
  pos <- try(grep(paste0("^", ws$trait[x], "$"), hmdb$synonyms, ignore.case = TRUE))
  
  eg1 <- try(hmdb[hmdb$name==hmdb$name[pos],])
  
  eg2 <- try(hmdb[grep(paste0("^", ws$trait[x], "$"), hmdb$name, ignore.case = TRUE),])
  
  eg3 <- rbind(eg1, eg2)
  
  eg <- eg3
  
  # eg <- eg3[!duplicated(eg3),]
  
  ws2 <- ws[ws$trait==ws$trait[x],]
  
  if(is(eg, "try-error")) {
    
    eg <- data.frame(secondary_accessions = NA, name = NA, synonyms = NA, uniprot = NA, exposure = ws2$exposure[1], trait = ws2$trait[1], Data = "DRAISMA_2015", stringsAsFactors = F)
    
  } else if(nrow(eg)!=0) {
    
    eg <- merge(eg, ws2, all = TRUE)
    eg$Data <- "DRAISMA_2015"
    
    eg <- data.frame(eg, stringsAsFactors = F)
    
  } else {
    
    eg <- data.frame(secondary_accessions = NA, name = NA, synonyms = NA, uniprot = NA, exposure = ws2$exposure[1], trait = ws2$trait[1], Data = "DRAISMA_2015", stringsAsFactors = F)
    
  }
  
  
}, mc.cores = parallel::detectCores(), mc.preschedule = FALSE)

# df2 <- do.call(rbind, df)

df2 <- bind_rows(df)

# df2$V1 <- NULL

df2$source <- with(df2, ifelse(is.na(uniprot), NA, "hmdb"))

# will need to remove norepinephrine link to Hexadecanoylcarnitine, not sure why this matched or assigned a synonyms in the xml file but not listed in the website > RAISE concern

# df3 <- df2[-grep("Error in", df2$exposure),]

# mets <- unique(df2$exposure) # 44 out of 163
# 
# unmapped_df <- ws[!ws$exposure %in% mets,] # manually mapping these

mets1 <- unique(df2$exposure[which(is.na(df2$source))])

mets2 <- unique(df2$exposure[which(df2$source=="hmdb")])

mets3 <- setdiff(mets1, mets2)

unmapped_df <- ws[ws$exposure %in% mets3,]

unmapped_df$uniprot <- ""

unmapped_df$uniprot <- as.list(unmapped_df$uniprot)

# manual mapping

mets <- c("carnitine", "Isoleucine", "Phosphatidylcholine", "Sphingomyeline", "Hexose", "Threonine", "Valine")

mets_go <- c("GO:0009437", "GO:0006549", "GO:0046470", "GO:0006684", "GO:0019318", "GO:0006566", "GO:0006573")

metsdf <- data.frame(metabolites = mets, goterms = mets_go, stringsAsFactors = F)

library(httr)
library(jsonlite)
library(xml2)

lapply(seq(metsdf$metabolites), function (x) {
  
  go <- gsub(":", "%3A", metsdf$goterms[x])
  
  print(go)
  
  requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/annotation/search?goId=", go, "&taxonId=9606&evidenceCode=ECO%3A0000269%2C%20ECO%3A0006056%2C%20ECO%3A0000250%2C%20ECO%3A0000245%2C%20ECO%3A0000303%2C%20ECO%3A0000305&limit=200")
  
  r <- GET(requestURL, accept("application/json"))
  
  # stop_for_status(r)
  
  json <- toJSON(content(r))
  
  eg <- fromJSON(json)
  
  eg2 <- eg$results
  
  eg3 <- as.character(unique(eg2$geneProductId))
  
  eg4 <- gsub('.*:(.*)','\\1', eg3)
  
  # metsdf$uniprot[x] <<- list(eg4)
  
  unmapped_df$uniprot[grep(metsdf$metabolites[x], unmapped_df$trait)] <<- list(eg4)
  
})

unmapped_df$secondary_accessions <- NA

unmapped_df$name <- NA

unmapped_df$synonyms <- NA

unmapped_df$source <- with(unmapped_df, ifelse(is.na(unmapped_df$uniprot), NA, "GO"))

unmapped_df$Data <- "DRAISMA_2015"

df3 <- unmapped_df %>% select(secondary_accessions, name, synonyms, uniprot, exposure, trait, Data, source)

df3$uniprot[which(lengths(df3$uniprot)==0)] <- NA_character_ # for threonine

# names(df3)[6] <- "trait"  

df <- rbind(df2, df3)  

saveRDS(df, "draisma_map_hmdb_go.rds")
