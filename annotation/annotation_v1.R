# Annotation script

# Unfiltered and partially annotated SNP dataset====

## 1. Load libraries====

library(dplyr)
library(furrr)
library(rols)
library(MendelianRandomization)
library(data.table)

## 2. Removing a faulty breast cancer GWAS (GCST007236), flipping IBD/Teslovich (2010) study MR betas & removing bxy=1 where exposure and outcome is the same ====

df <- readRDS("~/mr_results_snps/full_protein_MR_snp_v1_v2_feb2022.rds")

studies_to_remove <- c("GCST007799", "GCST007800", "GCST007236") # this step should be made redundant when GWAS catalog updates their corrected harmonised sumstats and are re-ingested by us.

df <- df[-which(df$outcome %in% studies_to_remove),]

studies_to_flip <- c("GCST004131", "GCST004132", "GCST004133", "GCST001725", "GCST001728", "GCST001729", "GCST000964", "GCST000758", "GCST000760", "GCST000755", "GCST000759")

df$bxy <- with(df, ifelse(outcome %in% studies_to_flip, -bxy, bxy))

df$bzy <- with(df, ifelse(outcome %in% studies_to_flip, -bzy, bzy))

df <- df[-which(df$bxy==1 & df$Data=="OLLI_2017"),] # removing self-Olli MR studies

##3. cis-trans annotations====

df$var_probe <- with(df, paste0(snp, ":", exposure)) ## getting the unique variant probe list

var_probe_unique <- data.frame(var_probe=unique(df$var_probe), stringsAsFactors = F)

var_probe_unique$exposure <- sapply(seq(var_probe_unique$var_probe), function (x) strsplit(var_probe_unique$var_probe[x], ":")[[1]][5])

## combining (and then loading) the full protein ensid map (this has ensid maps for proteins from 7 protein GWAS)

## v1 <- list.files("~/ensid_maps", full.names = T, pattern = "protein_ensid_map.rds", ignore.case = T) 

## map <- rbindlist(lapply(v1, readRDS))

## saveRDS(map, "~/ensid_maps/combined_protein_ensid_map_v1.rds")

map <- readRDS("~/ensid_maps/combined_protein_ensid_map_v1.rds")

df2 <- merge(var_probe_unique, map, by = "exposure")

df2$snp_chr <- sapply(seq(df2$var_probe), function (x) strsplit(df2$var_probe[x], ":")[[1]][1])

df2$snp_pos <- sapply(seq(df2$var_probe), function (x) strsplit(df2$var_probe[x], ":")[[1]][2])

df2$snp_pos <- as.numeric(df2$snp_pos)

df2$cis_trans <-  with(df2, ifelse((snp_chr==chr & snp_pos >= r1 & snp_pos <= r2), "cis", "trans"))

pos <- match(df$var_probe, df2$var_probe) # key to add new columns to the origina protein snp MR dataset

df$protein_trait <- df2$trait[pos]

df$ensid <- df2$ensid[pos]

df$cis_trans <- df2$cis_trans[pos]

df$var_probe <- NULL # don't need this anymore

rm("df2", "map", "pos", "var_probe_unique", "studies_to_flip") # clearing stuff

gc() # getting some memory back

df$exp_out_gsmr <- with(df, paste0(exposure, "_", outcome))

##4. Fetching colocs with gwas beta & perform single-SNP MRs where possible====

source("~/download_coloc_results_bq.R")

coloc$right_study <- gsub("SUN2018", "SUN_2018", coloc$right_study)

coloc$right_study <- gsub("OLLI_2016", "OLLI_2017", coloc$right_study)

coloc$right_study <- gsub("FOLKERSEN_2020", "FOLKERSEN_2017", coloc$right_study)

## merging coloc with main dataset using study IDs

df$key <- with(df, paste0(outcome, "_", Data, "_", ensid))

coloc$key <- with(coloc, paste0(left_study, "_", right_study, "_", right_gene_id))

df <- data.table(df) # converting to data.table results in faster merging

coloc <- data.table(coloc) # converting to data.table results in faster merging

df_coloc <- merge(df, coloc, by = "key", all = TRUE)

df_coloc$key <- NULL

studies_to_remove <- c("GCST007799", "GCST007800", "GCST007236") # this step should be made redundant when GWAS catalog updates their corrected harmonised sumstats and are re-ingested by us.

df_coloc <- df_coloc[-which(df_coloc$left_study %in% studies_to_remove),]

rm(list=setdiff(ls(), "df_coloc"))

gc()

## mapping coloc columns to the main data columns

### fixing gwas_beta columns (repositioning gwas_beta_from_or to gwas_beta columns)

df_coloc$gwas_beta <- with(df_coloc, ifelse(is.na(gwas_beta) & !is.na(gwas_beta_from_or), gwas_beta_from_or, gwas_beta))

df_coloc$gwas_se <- with(df_coloc, ifelse(is.na(gwas_se) & !is.na(gwas_se_from_orci), gwas_se_from_orci, gwas_se))

df_coloc[,c("gwas_beta_from_or", "gwas_se_from_orci")] <- NULL

df_coloc$coloc <- with(df_coloc, ifelse(!is.na(coloc_n_vars), "Yes", "No"))

df_coloc$exp_out_coloc <- with(df_coloc, ifelse(coloc=="Yes", paste0(right_study, "_", left_study, "_coloc"), NA))

### writing a function to map columns from coloc to main dataset

map_coloc <- function (df, cols) {
  
  # df is the name of the main dataset where the column values should be added
  
  # cols should provide variable names in the coloc column in the following order as strings: variant, beta_effect, beta_effect_se, beta_effect_pval, beta_outcome, beta_outcome_se, beta_outcome_pval, exposure_dataset, exposure_name, outcome_name, ensemble_id
  
  names_cols_main <- c("snp", "bzx", "bzx_se", "bzx_pval", "bzy", "bzy_se", "bzy_pval", "Data", "exposure", "outcome", "ensid") # add nsnp=1 and cis_trans="cis" later
  
  df <- get(df)
  
  # print(names(df))
  
  df2 <- data.frame(main_dataset_cols = names_cols_main, coloc_dataset_cols = cols, stringsAsFactors = F)
  
  # print(df2)
  
  sapply(seq(df2$main_dataset_cols), function (y) {
    
    var1 <- df2$main_dataset_cols[y]
    var2 <- df2$coloc_dataset_cols[y]
    
    # Mapping values
    
    pos1 <- which(is.na(df[,var1]) & !is.na(df[,var2]))
    
    df[,var1][pos1] <<- df[,var2][pos1]
    
    df[,"cis_trans"][pos1] <<- "cis"
    
    df[,"nsnp"][pos1] <<- 1 
    
    # df[,var2][pos1] <<- ""
    
  })
  
  df
  
}
  
df_coloc <- data.frame(df_coloc, stringsAsFactors = F)

df_coloc <- map_coloc("df_coloc", c("varid_left", "left_var_right_study_beta", "left_var_right_study_se", "left_var_right_study_pval", "gwas_beta", "gwas_se", "pval", "right_study", "right_phenotype", "left_study", "right_gene_id")) # column names should be in this order

gc()

### single-SNP MR with coloc var where possible

pos2 <- which(is.na(df_coloc[, "bxy"]) & !is.na(df_coloc[, "bzx"]) & !is.na(df_coloc[, "bzy"]))

eg <- df_coloc %>% select(bxy, bxy_se, bxy_pval, bzx, bzx_se, bzy, bzy_se) %>% slice(pos2)

sapply(seq(eg$bxy), function (x) {
  
  # pos <- pos2[x]
  
  bzx <- eg$bzx[x]
  bzx_se <- eg$bzx_se[x]
  bzy <- eg$bzy[x]
  bzy_se <- eg$bzy_se[x]
  
  # key <- df_coloc$key[pos]
  
  mr <- mr_ivw(mr_input(bx = bzx, bxse = bzx_se, by = bzy, byse = bzy_se))
  
  eg$bxy[x] <<- mr@Estimate
  eg$bxy_se[x] <<- mr@StdError
  eg$bxy_pval[x] <<- mr@Pvalue
  
  # mrdf <- data.frame(key, bxy, bxy_se, bxy_pval, stringsAsFactors = F)
  
})

#df_coloc$snp[pos2] <- eg$varid_left
df_coloc$bxy[pos2] <- eg$bxy
df_coloc$bxy_se[pos2] <- eg$bxy_se
df_coloc$bxy_pval[pos2] <- eg$bxy_pval

rm(list=setdiff(ls(), "df_coloc"))

gc()

##5. Assigning gene symbols to probe ensids====

genes <- unique(df_coloc$ensid)

source("~/ensid_hgnc.R")

pos <- match(df_coloc$ensid, hgnc$gene_id)

df_coloc$hgnc_protein <- hgnc$gene_name[pos]

rm(list=setdiff(ls(), "df_coloc"))

gc()

##6. Adding study EFOs====

### direct EFO and parent EFOs

source("~/downloading_studies_bq.R")

pos <- match(df_coloc$outcome, efo$study_id)

df_coloc$outcome_trait <- efo$trait_reported[pos]

df_coloc$outcome_trait_category <- efo$trait_category[pos]

df_coloc$outcome_trait_efo <- efo$trait_efos[pos]

df_coloc$n_initial <- efo$n_initial[pos]

df_coloc$n_cases <- efo$n_cases[pos]

rm(list=setdiff(ls(), "df_coloc"))

gc()

### parent efo

source("~/download_parent_efo_bq.R") # downloads parent efo

pefo2 <- pefo %>% group_by(study_id) %>% mutate(anc2 = list(paste(sort(unique(anc)),collapse=", "))) %>% ungroup() %>% select(study_id, anc2) %>% distinct()

pos <- match(df_coloc$outcome, pefo2$study_id)

df_coloc$parent_outcome_trait_efo <- pefo2$anc2[pos]

rm(list=setdiff(ls(), "df_coloc"))

gc()

# trms <- unique(unlist(df_coloc$outcome_trait_efo[which(nchar(df_coloc$outcome_trait_efo)!=0)]))

# source("~/download_parent_efo_bq.R") # downloads parent efo
# 
# pefo$parent_efo <- sapply(seq(pefo$id), function (x) pefo$anc[[x]]$element)

# trms <- trms[!which(nchar(trms)==0)]

# ont <- Ontology("efo")

# t <- term(ont, "EFO_0004574")

# trmslist <- parallel::mclapply(seq(trms), function (x) { # sometimes you will need to restart r session for this to work
#   
#   id <- trms[x]
#   
#   print(id)
#   
#   ont <- strsplit(id, "_")[[1]][1]
#   
#   print(ont)
#   
#   if((ont=="Orphanet")==TRUE) {
#     
#     ont <- "ORDO"
#     
#   }
#   
#   ontol <- Ontology(ont)
#   
#   trm <- term(ontol, id)
#   
#   if(is.null(parents(trm))==FALSE) {
#     
#     # p <- as(parents(trm), "data.frame")
#     
#     # parent_id <- list(p$id)
#     
#     # parent_label <- list(p$label)
#     
#     p <- parents(trm)
#     
#     p2 <- capture.output(p@x)
#     
#     parent_id <- gsub("\\$|\\`", "", p2[grep("\\$", p2)])
#     
#   } else {
#     
#     parent_id <- "no parent terms"
#     
#     # parent_label <- "no parent terms"
#     
#   }
#   
#   # df <- data.frame(id=id, pid=parent_id, plabel=parent_label, stringsAsFactors = F)
#   
#   df <- data.frame(id=id, pid=parent_id, stringsAsFactors = F)
#   
#   # names(df) <- c("id", "pid", "plabel")
#   
#   names(df) <- c("id", "pid")
#   
#   df
#   
# }, mc.cores = parallel::detectCores())
# 
# trmsdf <- do.call(rbind, trmslist)

# trmsdf2 <- trmsdf %>% tidyr::chop(cols = c("pid")) # can be unnested or unchopped later

# pefo2 <- pefo %>% group_by(id) %>% mutate(parent_efo2 = list(paste(sort(unique(parent_efo)),collapse=", "))) %>% ungroup() %>% select(id, parent_efo2)
# 
# pefo3 <- pefo2[-which(duplicated(pefo2)),]
# 
# eg <- df_coloc %>% select(outcome, outcome_trait_efo)
# 
# eg2 <- eg %>% tidyr::unchop(cols = c("outcome_trait_efo"), keep_empty = T)
# 
# eg3 <- eg[-which(dupl)]
# 
# pos <- match(eg2$outcome_trait_efo, pefo3$id)
# 
# eg2$parent_outcome_trait_efo <- pefo3$parent_efo2[pos]
# 
# eg3 <- eg2 %>% tidyr::chop(cols = c("outcome_trait_efo"))

# df_coloc2 <- df_coloc %>% tidyr::unchop(cols = c("outcome_trait_efo"), keep_empty = T)

# df_coloc3 <- df_coloc2[-which(duplicated(df_coloc2)),]

# efos <- unique(unlist(df_coloc$outcome_trait_efo[which(nchar(df_coloc$outcome_trait_efo)!=0)]))
# 
# pos1 <- match(efos, pefo3$id)
# 
# pefo4 <- data.frame(dfefos = pefo3$id[pos1], parent_efo = pefo3$parent_efo2[pos1], stringsAsFactors = F)

# pos <- df_coloc$outcome_trait_efo %in% list(pefo3$id)

# two-way matching strategy: 1) straightforward matching - works when there is only 1 outcome trait, 2) convoluted matching - works when there is more than 1 outcome trait that straightforward matching would result in NAs (I can do this using a single function but it takes a lot of time)

# quarters <- as.numeric(substr(names(sort(unlist(df_coloc$outcome_trait_efo))),1,1))
# 
# pos <- findInterval(match(b, unlist(a)), cumsum(c(0,sapply(a, length)))+1)
# 
# a <- df_coloc$outcome_trait_efo
# g <- rep(seq_along(a), sapply(a, length))
# pos <- g[match(unlist(a), pefo3$id),]
# aa <- unlist(a)
# au <- unique(aa)
# af <- factor(aa, levels=au)
# gg <- split(g, af)
# pos <- gg[match(au, pefo3$id)]


# pos <- match(df_coloc$outcome_trait_efo, pefo3$id)

# df_coloc$parent_outcome_trait_efo <- pefo3$parent_efo2[pos] # but does not match for multiple efos in list in df_coloc$outcome_trait_efo, doing this separately

# pos <- match(df_coloc$outcome_trait_efo, trmsdf2$id)

# df_coloc$parent_outcome_trait_efo <- trmsdf2$pid[pos] # but does not match for multiple efos in list in df_coloc$outcome_trait_efo, doing this separately

# pos_na <- which(is.na(df_coloc$parent_outcome_trait_efo))

# eg <- df_coloc[3000000,]
# 
# pos_null <- which(lengths(df_coloc$outcome_trait_efo) > 1)
# 
# parents_efo <- parallel::mclapply(seq(df_coloc$outcome_trait_efo)[pos_null], function (x) {
# 
#   # print(unlist(df_coloc$outcome_trait_efo[x]))
# 
#   r1 <- unlist(df_coloc$outcome_trait_efo[x])
# 
#   r2 <- pefo3$parent_efo2[match(r1, pefo3$id)]
# 
#   # print(unlist(r2))
# 
#   r3 <- unlist(r2)
# 
# }, mc.cores = parallel::detectCores())
# 
# df_coloc$parent_outcome_trait_efo[pos_na] <- parents_efo
# 
# class(df_coloc$parent_outcome_trait_efo) <- "list"

# saveRDS(df_coloc, "mr_prot_snp_unfiltered_dataset_v1_v2.rds")

##7. Get the fpred_max_label for all variants, assign PAVs to cis, v2g_trans for trans-variants that don't have an fpred_max_label====

snps <- gsub(":", "_", unique(df_coloc$snp))

source("~/SNP_vep_ld_annotation.R") # creates snpsdf (WILL NEED TO RE-WRITE THIS CODE)

snpsdf$varid_b38_index <- gsub("_", ":", snpsdf$varid_b38_index)

df_coloc$key <- with(df_coloc, paste0(snp, "_", ensid))

snpsdf$key1 <- with(snpsdf, paste0(varid_b38_index, "_", gene_id_index))

snpsdf$key2 <- with(snpsdf, paste0(varid_b38_index, "_", gene_id_tag))

posgene1 <- match(df_coloc$key, snpsdf$key1)

df_coloc$fpred_max_label_index <- snpsdf$fpred_max_label_index[posgene1]

posgene2 <- match(df_coloc$key, snpsdf$key2)

df_coloc$fpred_max_label_tag <- snpsdf$fpred_max_label_tag[posgene2]

df_coloc$pav_cis <- with(df_coloc, ifelse(cis_trans=="cis" & 
                              (!is.na(fpred_max_label_index) & 
                                 fpred_max_label_index %in% c("missense_variant", "stop_gained", "stop_lost", "start_lost", " transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant", "transcript_amplification", "inframe_insertion", "inframe_deletion", "protein_altering_variant", "TFBS_ablation")) | 
                              (!is.na(fpred_max_label_tag) &  fpred_max_label_tag %in% c("missense_variant", "stop_gained", "stop_lost", "start_lost", " transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant", "transcript_amplification", "inframe_insertion", "inframe_deletion", "protein_altering_variant", "TFBS_ablation")), "Yes", "No"))

df_coloc$key <- NULL

rm(list=setdiff(ls(), "df_coloc"))

gc()

##8. Get the v2g highest scored genes and apply them only to trans-pQTLs that don't have an fpred_max_label====

snps <- unique(df_coloc$snp)

source("~/SNP_v2g_trans_annotation.R")

pos <- match(df_coloc$snp, v2g$varid)

df_coloc$v2g <- v2g$gene_id[pos]

##9. Assigning gene symbols to v2g ensids====

genes <- unique(df_coloc$v2g)

source("~/ensid_hgnc.R")

pos <- match(df_coloc$v2g, hgnc$gene_id)

df_coloc$hgnc_v2g <- hgnc$gene_name[pos]

rm(list=setdiff(ls(), "df_coloc"))

gc()

##10. Saving unfiltered MR SNP file==== 

saveRDS(df_coloc, "~/mr_prot_snp_unfiltered_dataset_v1_v2.rds")

rm(list=setdiff(ls(), "df_coloc"))

gc()

# Filtered and Egger/WM-annotated SNP dataset ====

# df_coloc <- readRDS("~/mr_prot_snp_unfiltered_dataset_v1_v2.rds")

df_coloc <- df_coloc[which((!is.na(df_coloc$bxy_pval) & df_coloc$bxy_pval < 0.0005) | (is.na(df_coloc$bxy_pval) & !is.na(df_coloc$coloc_h4_h3) & df_coloc$coloc_h4_h3 > 1)),]

gc()

df_coloc$exp_out_gsmr_coloc <- with(df_coloc, ifelse(!is.na(exp_out_gsmr) & !is.na(exp_out_coloc), paste0(exp_out_gsmr, "_", exp_out_coloc, "_", ensid, "_", varid_left), ifelse(is.na(exp_out_gsmr) & !is.na(exp_out_coloc), paste0(exp_out_coloc, "_", ensid, "_", varid_left), exp_out_gsmr)))


##1. MR-Egger and MR-weighted median====

# df_coloc$exp_out <- with(df_coloc, paste0(exposure, "_", outcome))
# 
# df_coloc$exp_out[which(df_coloc$coloc=="Yes")] <- paste0(df_coloc$exp_out[which(df_coloc$coloc=="Yes")] , "_coloc")

mr_egger_wm <- function (x) { # this takes the most time, ~ 2-3 hours
  
  if(x[, "nsnp"][1][[1]] >= 3) {
    
    eg <- with(x, mr_allmethods(mr_input(bx = bzx, bxse = bzx_se, by = bzy, byse = bzy_se), method = "main"))
    
    # eg <- with(x, mr_allmethods(mr_input(bx = bx1, bxse = bxse1, by = by1, byse = byse1), method = "main"))
    
    mr_egger <- eg@Values$Estimate[eg@Values$Method=="MR-Egger"]
    
    mr_egger_se <- eg@Values$`Std Error`[eg@Values$Method=="MR-Egger"]
    
    mr_egger_p <- eg@Values$`P-value`[eg@Values$Method=="MR-Egger"]
    
    mr_wm <- eg@Values$Estimate[eg@Values$Method=="Weighted median"]
    
    mr_wm_se <- eg@Values$`Std Error`[eg@Values$Method=="Weighted median"]
    
    mr_wm_p <- eg@Values$`P-value`[eg@Values$Method=="Weighted median"]
    
  } else {
    
    mr_egger <- NA
    
    mr_egger_se <- NA
    
    mr_egger_p <- NA
    
    mr_wm <- NA
    
    mr_wm_se <- NA
    
    mr_wm_p <- NA
    
  }
  
  out <- data.frame(mr_egger, mr_egger_se, mr_egger_p, mr_wm, mr_wm_se, mr_wm_p)
  
  out2 <- cbind.data.frame(x, out)
  
  return(out2)
  
}

# library(furrr)
future::plan(multisession)

df <- df_coloc %>% split(.$exp_out_gsmr_coloc) %>% furrr::future_map_dfr(mr_egger_wm, .progress = TRUE)

rm(list=setdiff(ls(), "df"))

gc()

##2. Saving filtered MR SNP file==== 

saveRDS(df, "~/mr_prot_snp_filtered_dataset_v1_v2.rds")

# Filtered and annotated MR dataset====

##1. new MR columns====

df <- readRDS("~/mr_prot_snp_filtered_dataset_v1_v2.rds")

# df$exp_out_gsmr_coloc <- with(df, ifelse(!is.na(exp_out_gsmr) & !is.na(exp_out_coloc), paste0(exp_out_gsmr, "_", exp_out_coloc), ifelse(is.na(exp_out_gsmr) & !is.na(exp_out_coloc), exp_out_coloc, exp_out_gsmr)))


dfmr <- df %>% group_by(exp_out_gsmr_coloc) %>% 
  mutate(v2g_mr = list(paste(sort(unique(v2g)),collapse=", ")), 
         hgnc_v2g_mr = list(paste(sort(unique(hgnc_v2g)),collapse=", ")), 
         cis_trans_mr = ifelse((all(cis_trans=="cis") | all(cis_trans=="trans")), cis_trans, "mixed"), 
         pav_cismr = ifelse(any(!is.na(pav_cis)) & any(pav_cis=="Yes") & all(cis_trans=="cis"), "Yes", NA),
         fpred_max_label_index_mr = list(paste(sort(unique(fpred_max_label_index)),collapse=", ")),
         fpred_max_label_tag_mr = list(paste(sort(unique(fpred_max_label_tag)),collapse=", "))) %>%   ungroup() %>% 
  select(Data, exposure, protein_trait, hgnc_protein, ensid, outcome, outcome_trait, exp_out_gsmr, exp_out_coloc, exp_out_gsmr_coloc, varid_left, nsnp, n_cases, n_initial, outcome_trait_efo, parent_outcome_trait_efo, cis_trans_mr, coloc, coloc_n_vars, coloc_h1, coloc_h2, coloc_h3, coloc_h4, coloc_h4_h3, fpred_max_label_index_mr, fpred_max_label_tag_mr, bxy, bxy_se, bxy_pval, mr_egger, mr_egger_se, mr_egger_p, mr_wm, mr_wm_se, mr_wm_p, pav_cismr, v2g_mr, hgnc_v2g_mr) %>% 
  distinct()

# Corrections 1: after grouping by exp_out of MR filtered dataset, removing duplicate rows that have NULL parent EFO

# dfmr <- readRDS("~/mr_prot_filtered_dataset_v1_v2.rds")

# expout <- unique(dfmr$exp_out_gsmr_coloc)

# dflst <- lapply(seq(expout), function (x) { # remving duplicate rows with null parent EFOs
#   
#   # print(expout[x])
#   
#   df <- dfmr[which(dfmr$exp_out_gsmr_coloc==expout[x]),]
#   
#   dup <- any(duplicated(df$exp_out_gsmr_coloc))
#   
#   if (dup==TRUE) {
#   
#   names(df$parent_outcome_trait_efo) <- seq_along(df$parent_outcome_trait_efo)
#   
#   # df2 <- df[lapply(df$parent_outcome_trait_efo, length) > 0,]
#   
#   df <- df[!sapply(df$parent_outcome_trait_efo, is.null),]
#   
#   } else {
#     
#     df
#     
#   }
#   
# })
# 
# dfmr2 <- data.table::rbindlist(dflst, fill = TRUE)

##2. separating p-value mantissa from its exponent====

mantissa_exponent_zero <- function (df, col) { # a better function than above - creates both mantissa and exponent + deals with 0 p-values
  
  # print(length(col))
  
  for (i in 1:length(col)) {
    
    # print(col[i])
    
    beta <- gsub("_pval|_p", "", col[i])
    
    se <- gsub("_pval|_p", "_se", col[i])
    
    # df <- get(df)
    
    # print(names(df))
    
    lst <- lapply(seq(df[,col[i]]), function (x) {
      
      if(df[,col[i]][x]==0 & !is.na(df[,col[i]][x])) {
        
        eg <- 2*Rmpfr::pnorm(Rmpfr::mpfr(abs(df[,beta][x]/df[,se][x]), precBits=100), lower.tail=FALSE, log.p = FALSE)
        
        mantissa <- as.numeric(strsplit(gsub("\\[1\\] ", "", as.character(capture.output(eg))[2]), "e")[[1]][1])
        
        exponent <- as.numeric(strsplit(gsub("\\[1\\] ", "", as.character(capture.output(eg))[2]), "e")[[1]][2])
        
        lst <- c(mantissa, exponent)
        
      } else {
        
        # col <- paste0(col[i], "_mantissa")
        
        # mantissa <- df[,col][x]
        
        mantissa <- as.numeric(strsplit(format(df[,col[i]][x], scientific = TRUE, digits = 5), "e")[[1]][1])
        
        exponent <- as.numeric(strsplit(format(df[,col[i]][x], scientific = TRUE, digits = 5), "e")[[1]][2])
        
        lst <- c(mantissa, exponent)
        
      }
      
    })
    
    lst <- unlist(lst)
    
    # print(col[i])
    
    # print(lst[c(TRUE, FALSE)])
    
    mantissa <- lst[c(TRUE, FALSE)] # odd rows
    
    exponent <- lst[!c(TRUE, FALSE)] # even rows
    
    varname1 <- paste0(col[i], "_mantissa")
    
    varname2 <- paste0(col[i], "_exponent")
    
    df[[varname1]] <- mantissa
    
    df[[varname2]] <- exponent
    
    df
    
  }
  
  df
  
}

dfmr <- data.frame(dfmr, stringsAsFactors = F)

dfmr <- mantissa_exponent_zero(dfmr, c("bxy_pval", "mr_egger_p", "mr_wm_p"))

##3. Getting odds ratio and CIs====

orci <- function (df, col) { 
  
  for (i in 1:length(col)) {
    
    beta <- col[i]
    
    se <- paste0(col[i], "_se")
    
    lst <- lapply(seq(df[,col[i]]), function (x) {
      
      if(!is.na(df[,col[i]][x])) {
        
        # obtain 95% CI from beta
        
        beta <- df[,beta][x]
        se <- df[,se][x]
        
        beta_lci <- beta-1.96*se
        beta_uci <- beta+1.96*se
        
        # OR
        
        or <- exp(beta)
        
        # 95% CI of OR
        
        or_lci <- exp(beta_lci)
        
        or_uci <- exp(beta_uci)
        
        lst <- c(or, or_lci, or_uci, beta_lci, beta_uci)
        
        
      } else {
        
        or <- NA
        or_lci <- NA
        or_uci <- NA
        beta_lci <- NA 
        beta_uci <- NA
        
        lst <- c(or, or_lci, or_uci, beta_lci, beta_uci)
        
      }
      
    })
    
    lst <- unlist(lst)
    
    or <- lst[c(TRUE, FALSE, FALSE, FALSE, FALSE)]
    
    or_lci <- lst[c(FALSE, TRUE, FALSE, FALSE, FALSE)]
    
    or_uci <- lst[c(FALSE, FALSE, TRUE, FALSE, FALSE)]
    
    beta_lci <- lst[c(FALSE, FALSE, FALSE, TRUE, FALSE)]
    
    beta_uci <- lst[c(FALSE, FALSE, FALSE, FALSE, TRUE)]
    
    varname1 <- paste0("or_", col[i])
    
    varname2 <- paste0("or_lci_", col[i])
    
    varname3 <- paste0("or_uci_", col[i])
    
    varname4 <- paste0(col[i], "_lci")
    
    varname5 <- paste0(col[i], "_uci")
    
    df[[varname1]] <- or
    
    df[[varname2]] <- or_lci
    
    df[[varname3]] <- or_uci
    
    df[[varname4]] <- beta_lci
    
    df[[varname5]] <- beta_uci
    
  }
  
  df
  
}

dfmr <- orci(dfmr, c("bxy", "mr_egger", "mr_wm"))

removing_or_for_qtraits <- function (x, data) {
  
  with(data, ifelse(is.na(n_cases), NA, x))
  
}

dfmr[,grep("or_", names(dfmr), value = T)] <- lapply(dfmr[,grep("or_", names(dfmr), value = T)], removing_or_for_qtraits, data=dfmr)

##4. Saving filtered MR file

saveRDS(dfmr, "~/mr_prot_filtered_dataset_v1_v2.rds")

rm(list=ls())

gc()

# Exporting datasets to cloud (optional)====

df_snp_unfiltered <- readRDS("~/mr_prot_snp_unfiltered_dataset_v1_v2.rds")
df_snp_filtered <- readRDS("~/mr_prot_snp_filtered_dataset_v1_v2.rds")
df_mr_filtered <- readRDS("~/mr_prot_filtered_dataset_v1_v2.rds")

all <- ls(pattern = "df")

# convert list columns to character (because we are saving as tsv file)

df <- lapply(seq(all), function (y) {
  
  require(dplyr)
  
  print(all[y])
  
  df <- get(all[y])
  
  df2 <- df %>% rowwise() %>% mutate(across(where(is.list), toString))
  
  df2
  
})

names(df) <- all

rm(list=setdiff(ls(), "df"))

gc()

list2env(df, globalenv())

rm(df)

# save as tsv.gz and exporting to google cloud

version <- "v2"

all <- ls(pattern = "df")

sapply(seq(all), function (y) {
  
  df <- get(all[y])
  
  all2 <- gsub("df", "prot", all[y])
  
  dfn <- all2
  
  dfn_zipped <- paste0(dfn, ".tsv.gz")
  
  write.table(df, file=gzfile(paste0("~/", dfn_zipped)), sep = "\t", quote = F, row.names = F)
  
  print(paste0(dfn_zipped, " is created"))
  
  system(paste0("gsutil cp -r ", dfn_zipped, " gs://genetics-portal-dev-mr/mr_results_", version, "/"))
  
  print(paste0(dfn_zipped, " has been exported to GC"))
  
  
})

# save as just tsv

all2 <- all[-3] # removing unfiltered

sapply(seq(all), function (y) {
  
  df <- get(all2[y])
  
  all3 <- gsub("df", "prot", all2[y])
  
  dfn <- all3
  
  dfn <- paste0(dfn, ".tsv")
  
  write.table(df, file=paste0("~/", dfn), sep = "\t", quote = F, row.names = F)
  
  print(paste0(dfn, " is created"))
  
  system(paste0("gsutil cp -r ", dfn, " gs://genetics-portal-dev-mr/mr_results_", version, "/"))
  
  print(paste0(dfn, " has been exported to GC"))
  
})



# Corrections 2: remove GCST007799, GCST007800 (From both MR and coloc), and GCST007236 (from coloc) - this should be made redundant when re-ingested corrected harmonised sumstats:

# studies_to_remove <- c("GCST007799", "GCST007800", "GCST007236")
# 
# # unfiltered
# 
# df1 <- readRDS("~/mr_prot_snp_unfiltered_dataset_v1_v2.rds")
# 
# df2 <- df1[-which(df1$outcome %in% studies_to_remove),]
# 
# saveRDS(df2, "~/mr_prot_snp_unfiltered_dataset_v1_v2.rds")
# 
# rm(df1, df2)
# 
# gc()
# 
# # filtered
# 
# dfmr <- readRDS("~/mr_prot_filtered_dataset_v1_v2.rds")
# dfsnp <- readRDS("~/mr_prot_snp_filtered_dataset_v1_v2.rds")
# 
# dfmr <- dfmr[-which(dfmr$outcome %in% studies_to_remove),]
# dfsnp <- dfsnp[-which(dfsnp$outcome %in% studies_to_remove),]
# 
# saveRDS(dfmr, "~/mr_prot_filtered_dataset_v1_v2.rds")
# saveRDS(dfsnp, "~/mr_prot_snp_filtered_dataset_v1_v2.rds")


# eg <- dfmr2[which(dfmr2$coloc=="Yes" & dfmr2$nsnp!=1),]
# 
# eg$exp_out2 <- paste(eg$exp_out, "_", eg$varid_left)
# 
# eg$exp_out <- NULL


# Saving new files

# dfsnp <- readRDS("~/mr_prot_snp_filtered_dataset_v1_v2.rds")
# 
# dfsnp2 <- dfsnp[which(dfsnp$exp_out %in% expout2),]
# 
# saveRDS(dfmr2, "~/mr_prot_filtered_dataset_v1_v2.rds")
# 
# saveRDS(dfsnp2, "~/mr_prot_snp_filtered_dataset_v1_v2.rds")


# To deploy to shiny app====

# system("mkdir mr_app2/")
# system("cp mr_prot_filtered_dataset_v1_v2.rds mr_app2/")
# system("cp mr_prot_snp_filtered_dataset_v1_v2.rds mr_app2/")


# change in server.R df <- readRDS(<MR file>)
# edit ui.r appropriately

# rsconnect::deployApp("~/mr_app") # use when ready


