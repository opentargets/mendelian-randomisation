# extracting unfiltered MR data from unfiltered SNP data for analysis

library(dplyr)

df <- readRDS("~/mr_prot_snp_unfiltered_dataset_v1_v2.rds")

df$exp_out_gsmr_coloc <- with(df, ifelse(!is.na(exp_out_gsmr) & !is.na(exp_out_coloc), paste0(exp_out_gsmr, "_", exp_out_coloc, "_", ensid, "_", varid_left), ifelse(is.na(exp_out_gsmr) & !is.na(exp_out_coloc), paste0(exp_out_coloc, "_", ensid, "_", varid_left), exp_out_gsmr)))

# df$exp_out_gsmr_coloc <- with(df, ifelse(!is.na(exp_out_gsmr) & !is.na(exp_out_coloc), paste0(exp_out_gsmr, "_", exp_out_coloc), ifelse(is.na(exp_out_gsmr) & !is.na(exp_out_coloc), exp_out_coloc, exp_out_gsmr)))


dfmr <- df %>% group_by(exp_out_gsmr_coloc) %>% 
  mutate(cis_trans_mr = ifelse((all(cis_trans=="cis") | all(cis_trans=="trans")), cis_trans, "mixed"), 
         pav_cismr = ifelse(any(!is.na(pav_cis)) & any(pav_cis=="Yes") & all(cis_trans=="cis"), "Yes", NA),
         fpred_max_label_index_mr = list(paste(sort(unique(fpred_max_label_index)),collapse=", ")),
         fpred_max_label_tag_mr = list(paste(sort(unique(fpred_max_label_tag)),collapse=", "))) %>%   ungroup() %>% 
  select(Data, exposure, protein_trait, hgnc_protein, ensid, outcome, outcome_trait, exp_out_gsmr, exp_out_coloc, exp_out_gsmr_coloc, varid_left, nsnp, n_cases, n_initial, outcome_trait_efo, parent_outcome_trait_efo, cis_trans_mr, coloc, coloc_h4, coloc_h4_h3, fpred_max_label_index_mr, fpred_max_label_tag_mr, bxy, bxy_se, bxy_pval, pav_cismr) %>% 
  distinct()

saveRDS(dfmr, "~/mr_prot_unfiltered_dataset_v1_v2_without_egger.rds")
