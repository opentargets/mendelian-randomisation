# aligning metabolomic maps to protein maps (cols: exposure, trait, ensid, chr, tss, r1, r2)

# targets df

f <- dir("~/targets", full.names = T)

f <- f[-grep("_SUCCESS", f)]

lines <- lapply(seq(f), function (x) {
  
  l <- readLines(f[x])
  
  m <- jsonlite::stream_in(textConnection(l), verbose = F)
  
})

targets <- dplyr::bind_rows(lines)

targets2 <- targets[-which(is.na(targets$proteinAnnotations$id)),]

targets3 <- data.frame(ensid = targets2$id, uniprot = targets2$proteinAnnotations$id, stringsAsFactors = F)

# genes/tss/chr

library(bigrquery)

project <- "open-targets-genetics" # replace this with your project ID 

qry1 <- "SELECT
  *
FROM
  `210608.genes`"

genes <- bq_project_query(project, qry1) %>% bq_table_download()

# shin map

shin <- readRDS("shin_map_hmdb.rds")

pos <- match(shin$uniprot, targets2$proteinAnnotations$id)

shin$ensid <- targets2$id[pos]

shin2 <- shin %>% select(exposure, trait, uniprot, ensid, source)

shin3 <- shin2[which(!duplicated(shin2)),]

# group by exposure, and ifelse, if any row contains an ensid, retain the row(s), but if else no row contains an ensid, retain only the first (or any) row

shin4 <- shin3 %>% group_by(exposure) %>% mutate(delete1 = ifelse(is.na(ensid) & is.na(source), "yes", "no"), delete2 = ifelse((length(ensid)==1 & delete1=="yes"), "no", "yes"))

shin5 <- shin4[-which(shin4$delete1=="yes" & shin4$delete2=="yes"),]

shin5[,c("delete1", "delete2")] <- NULL

# adding tss and r1/r2 parameters

pos2 <- match(shin5$ensid, genes$gene_id)

shin5$chr <- genes$chr[pos2]

shin5$tss <- genes$tss[pos2]

shin5$r1 <- shin5$tss-1000000

shin5$r2 <- shin5$tss+1000000

# shin4 <- shin5[-which(is.na(shin5$chr)),]

# SHIN exposure format > formatted_<exposure>_metal_pos_hm_qc

# shin4$exposure <- paste0("formatted_", shin4$exposure, "_metal_pos_hm_qc")

shin5$exposure <- paste0("formatted_", shin5$exposure, "_metal_pos_hm_qc")

saveRDS(shin5, "~/ensid_maps/shin_metabolite_ensid_map.rds")

# kettunen map

kettunen <- readRDS("kettunen_map_hmdb_go.rds")

kettunen_all <- tidyr::unnest(kettunen, cols = c("uniprot"))

pos <- match(kettunen_all$uniprot, targets2$proteinAnnotations$id)

kettunen_all$ensid <- targets2$id[pos]

kettunen2 <- kettunen_all %>% select(exposure, trait, uniprot, ensid, source)

kettunen3 <- kettunen2[!duplicated(kettunen2),]

# group by exposure, and ifelse, if any row contains an ensid, retain the row(s), but if else no row contains an ensid, retain only the first (or any) row

kettunen4 <- kettunen3 %>% group_by(exposure) %>% mutate(delete1 = ifelse(is.na(ensid) & is.na(source), "yes", "no"), delete2 = ifelse((length(ensid)==1 & delete1=="yes"), "no", "yes"))

kettunen5 <- kettunen4[-which(kettunen4$delete1=="yes" & kettunen4$delete2=="yes"),]

kettunen5[,c("delete1", "delete2")] <- NULL

# adding tss and r1/r2 parameters

pos2 <- match(kettunen5$ensid, genes$gene_id)

kettunen5$chr <- genes$chr[pos2]

kettunen5$tss <- genes$tss[pos2]

kettunen5$r1 <- kettunen5$tss-1000000

kettunen5$r2 <- kettunen5$tss+1000000

# kettunen4 <- kettunen5[-which(is.na(kettunen5$chr)),]

# KETTUNEN exposure format > formatted_Summary_statistics_MAGNETIC_<exposure>_hm_qc

# kettunen4$exposure <- paste0("formatted_Summary_statistics_MAGNETIC_", kettunen4$exposure, "_hm_qc")

kettunen5$exposure <- paste0("formatted_Summary_statistics_MAGNETIC_", kettunen5$exposure, "_hm_qc")

kettunen5$exposure <- gsub("\\.", "_", kettunen5$exposure) # an additional step for kettunen

saveRDS(kettunen5, "~/ensid_maps/kettunen_metabolite_ensid_map.rds")

# draisma

draisma <- readRDS("draisma_map_hmdb_go.rds")

draisma_all <- tidyr::unnest(draisma, cols = c("uniprot"))

pos <- match(draisma_all$uniprot, targets2$proteinAnnotations$id)

draisma_all$ensid <- targets2$id[pos]

draisma2 <- draisma_all %>% select(exposure, trait, uniprot, ensid, source)

draisma3 <- draisma2[!duplicated(draisma2),]

# group by exposure, and ifelse, if any row contains an ensid, retain the row(s), but if else no row contains an ensid, retain only the first (or any) row

draisma4 <- draisma3 %>% group_by(exposure) %>% mutate(delete1 = ifelse(is.na(ensid) & is.na(source), "yes", "no"), delete2 = ifelse((length(ensid)==1 & delete1=="yes"), "no", "yes"))

draisma5 <- draisma4[-which(draisma4$delete1=="yes" & draisma4$delete2=="yes"),]

draisma5[,c("delete1", "delete2")] <- NULL

# adding tss and r1/r2 parameters

pos2 <- match(draisma5$ensid, genes$gene_id)

draisma5$chr <- genes$chr[pos2]

draisma5$tss <- genes$tss[pos2]

draisma5$r1 <- draisma5$tss-1000000

draisma5$r2 <- draisma5$tss+1000000

# draisma4 <- draisma5[-which(is.na(draisma5$chr)),]

# DRAISMA exposure format > formatted_MetaAnalysis_<exposure>_HetISqLt75_HetDfGt0_HWEPValGt1EMinus6_hm_qc

# draisma4$exposure <- paste0("formatted_MetaAnalysis_", draisma4$exposure, "_HetISqLt75_HetDfGt0_HWEPValGt1EMinus6_hm_qc")
# 
# draisma4$exposure <- gsub(":| |-", "_", draisma4$exposure)
# 
# draisma4$exposure <- gsub("\\(|\\)", "", draisma4$exposure)

draisma5$exposure <- paste0("formatted_MetaAnalysis_", draisma5$exposure, "_HetISqLt75_HetDfGt0_HWEPValGt1EMinus6_hm_qc")

draisma5$exposure <- gsub(":| |-", "_", draisma5$exposure)

draisma5$exposure <- gsub("\\(|\\)", "", draisma5$exposure)

saveRDS(draisma5, "~/ensid_maps/draisma_metabolite_ensid_map.rds")
