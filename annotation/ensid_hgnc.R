# Using BigQuery (genomics-public-data:linkage_disequilibrium_1000G_phase_3.super_pop_EUR & open-targets-genetics:210608.variants)  for LD matrix construction

# updated version feb2022: not using v2g to assign the top gene - using vep instead

library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID

qry1 <- paste(c("SELECT 
    gene_id,
    gene_name
FROM 
    `open-targets-genetics.210608.genes`
WHERE 
    gene_id IN (", toString(shQuote(genes)) ,")"), collapse = "")

hgnc <- bq_project_query(project, qry1) %>% bq_table_download(page_size = 30000) # toggle page_size if table is too large to be parsed

# df4$r2 <- df4$corr^2

# pav <- c("missense_variant", "stop_gained", "stop_lost", "start_lost", " transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant", "transcript_amplification", "inframe_insertion", "inframe_deletion", "protein_altering_variant", "TFBS_ablation") # https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html

# df4$ld_pav <- with(df4, ifelse(most_severe_consequence %in% pav & r2 >= 0.6, "Yes", "No")) # NEED TO FIX THIS

# df4 <- df4[which(!duplicated(df4$varid_b38_index)),]

# snpsdf <- df4 %>% select(varid_b38, rs_id, most_severe_consequence, gene_id_topv2g, ld_pav)
