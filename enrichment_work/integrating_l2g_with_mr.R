# annotate filtered MR table with L2G data

## Load filtered MR data

df_mr_filtered <- data.table::fread("~/prot_mr_filtered.tsv", stringsAsFactors = F)

df_mr_filtered$outcome_trait_efo[df_mr_filtered$outcome_trait_efo==""] <- NA

df_mr <- df_mr_filtered %>% tidyr::separate_rows(outcome_trait_efo, sep = ",") 

df_mr$key <- with(df_mr, paste0(ensid, "_", outcome))

## Load L2G table

project <- "open-targets-genetics" # replace this with your project ID

qry <- "SELECT DISTINCT
study_id,
CONCAT(chrom, ':', pos, ':', ref, ':', alt) as snp_l2g,
gene_id,
y_proba_full_model as l2g_score
FROM 
`open-targets-genetics-dev.genetics_dev.locus2gene`" # NEED TO PRE-JOIN ON EFO FIRST

l2g_tbl <- bq_project_query(project, qry) %>% bq_table_download()

l2g_tbl$key <- with(l2g_tbl, paste0(gene_id, "_", study_id))

l2g_tbl2 <- l2g_tbl %>% group_by(key) %>% filter(l2g_score==max(l2g_score)) %>% select(key, snp_l2g, l2g_score)

# merging MR with L2G

mr_l2g <- merge(df_mr, l2g_tbl2, by = "key", all.x = T)

# what proportion of MR have L2G < 0.5?

mr_no_l2g <- mr_l2g %>% filter(is.na(l2g_score))

mr_l2g_g0.5 <- mr_l2g %>% filter(l2g_score > 0.5)

mr_l2g_l0.5 <- mr_l2g %>% filter(l2g_score < 0.5)

length(unique(mr_l2g_l0.5$key)) # MR can help prioritise 1056 targets for 464 EFO traits where the L2G < 0.5. This set represents 9080 unique target-trait pairs

mr_l2g_l0.5_cases <- mr_l2g_l0.5[which(!is.na(mr_l2g_l0.5$n_cases)),]

mr_conf <- mr_l2g_l0.5_cases[which(mr_l2g_l0.5_cases$coloc_h4 > 0.8),]

# will need to do above with EFO