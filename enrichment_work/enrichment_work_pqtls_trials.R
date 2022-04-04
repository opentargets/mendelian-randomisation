# Can pQTL MR-coloc associations predict of success and failure of clinical trials?

# Approach: 
#1. use unique target-trait pairs to match between pQTL MR-coloc and drug datasets
#2. trait = parent EFOs
#3. test = fishers exact

#Q. This appraoch via parent EFOs duplicates the same studies due to multiple parent EFOs assigned to the same study (for enrichment analysis, this shouldn't be a problem as the duplications will affect counts in all cells of the 2x2 table)

library(dplyr)
library(forestplot)
library(ggplot2)

## MR dataset====

df <- readRDS("~/mr_prot_unfiltered_dataset_v1_v2_without_egger.rds")

# df$outcome_trait_efo <- NULL # will be using parent EFOs, so don't need this for now

df2 <- df %>% tidyr::unchop(cols = c("parent_outcome_trait_efo")) %>% tidyr::separate_rows("parent_outcome_trait_efo", sep = ", ")

# df3 <- df2[-which(df2$parent_outcome_trait_efo %in% c("EFO_0000651", "HP_0000118")),] # removing "phenotype" labels as this is non-specific

df3 <- df2[-which(df2$outcome_trait_category %in% c("phenotype", "Uncategorised", "injury, poisoning, or other complication")),] # removing non-specific labels

df3$key <- with(df3, paste0(ensid, "_", parent_outcome_trait_efo)) # key to merge with drug dataset

# df3$key <- with(df3, paste0(ensid, "_", outcome_trait_efo))

df4 <- df3 %>% group_by(key) %>% filter(is.na(bxy_pval) | bxy_pval == min(bxy_pval)) # within group, there will be multiple target-trait associations reported by multiple protein or outcome studies, so selecting the one with the lowest p-val (also retaining NA pvals as they may be a valid coloc association for which we were not able to perform MR analysis most likely due to missing exposure or outcome estimates)

df5 <- df4[-which(duplicated(df4$key)),] # retaining only unique target-trait pairs

df5$mrcoloc_association <- with(df5, ifelse(bxy_pval < 0.0005 | (bxy_pval < 0.0005 & coloc_h4_h3 > 1), 1, 0))

df5 <- data.table::data.table(df5) # converting to data.table as this can help is merge faster

rm(df, df2, df3, df4)

gc()

## Successful trials enrichment====

source("~/download_notstopped_trials_bq.R")

drugs$key <- with(drugs, paste0(ensid, "_", parent_efo))

# drugs$diseaseFromSourceMappedId <- NULL

# drugs <- drugs[-which(duplicated(drugs)),] 

# new indicator variables

drugs$phase_rev <- gsub("[^0-9.-]", "", drugs$Phase)

drugs$phase_rev <- with(drugs, ifelse(phase_rev=="12", "2", phase_rev)) # if phase 1/phase 2, label as phase 2

drugs$phase_rev <- with(drugs, ifelse(phase_rev=="23", "3", phase_rev)) # if phase 2/phase 3, label as phase 3

drugs$phase_rev <- as.integer(drugs$phase_rev)

# without a target-trait group, selecting the trials reporting success as the marker of success

drugs2 <- drugs %>% group_by(key) %>% filter(phase_rev==max(phase_rev)) %>% mutate(success = ifelse(phase_rev==4, 1, 0)) %>% select(key, trait_category, diseaseFromSourceMappedId, Phase, CT, Status, chembl_id, success)

drugs3 <- drugs2[-which(duplicated(drugs2$key)),]

rm(drugs, drugs2)

drugs3 <- data.table::data.table(drugs3) # converting to data.table as this can help is merge faster

# merge dataset by ensid and parent efo represented by the key

df_drugs <- merge(df5, drugs3, by = "key", all = TRUE)

rm(drugs3)

# indicator variables for enrichment analysis

df_drugs$genetic_association <- with(df_drugs, ifelse(mrcoloc_association==1, "1_Evidence of MR-coloc association", "2_No evidence of MR-coloc association"))

df_drugs$trial_status <- with(df_drugs, ifelse(success==1, "1_Phase 4", "2_Not Phase 4"))

# bar chart

eg_success <- df_drugs[with(df_drugs, which(mrcoloc_association==1 & success==1)),]

success_cat <- ggplot(eg_success, aes(x=reorder(trait_category, trait_category,
                         function(x)-length(x)))) + geom_bar() + theme(axis.text.x = element_text(angle = 45, vjust = 0.6))

# 2 x 2 table

success <- with(df_drugs, table(genetic_association, trial_status))

# fishers exact test

enrichment_success <- success %>% fisher.test()

## Terminated trials enrichment====

source("~/download_stopped_trials_bq.R")

drugs$key <- with(drugs, paste0(ensid, "_", parent_efo))

# drugs$key <- with(drugs, paste0(ensid, "_", diseaseFromSourceMappedId))

# drugs$diseaseFromSourceMappedId <- NULL
# 
# drugs$diseaseFromSource <- NULL
# 
# drugs$trait_category <- NULL
# 
# drugs <- drugs[-which(duplicated(drugs)),] 

# without a target-trait group, selecting the trials reporting success as the marker of success

drugs2 <- drugs %>% group_by(key) %>% mutate(failure = ifelse(Reason=="Safety_Sideeffects", 1, 0))

drugs3 <- drugs2[-which(duplicated(drugs2$key)),]

rm(drugs, drugs2)

drugs3 <- data.table::data.table(drugs3) # converting to data.table as this can help is merge faster

# merge dataset by ensid and parent efo represented by the key

df_drugs <- merge(df5, drugs3, by = "key", all = TRUE)

rm(drugs3)

# indicator variables for enrichment analysis

df_drugs$genetic_association <- with(df_drugs, ifelse(mrcoloc_association==1, "1_Evidence of MR-coloc association", "2_No evidence of MR-coloc association"))

df_drugs$trial_status <- with(df_drugs, ifelse(failure==1, "1_Terminated due to safety issues", "2_No evidence of termination due to safety issues"))

# bar chart

eg_failure <- df_drugs[with(df_drugs, which(mrcoloc_association==1 & failure==1)),]

# failure_cat <- ggplot(eg_failure, aes(x=reorder(trait_category, trait_category,
#                                         function(x)-length(x)))) + geom_bar() + theme(axis.text.x = element_text(angle = 45, vjust = 0.6))

# 2 x 2 table

failure <- with(df_drugs, table(genetic_association, trial_status))

# fishers exact test

enrichment_failure <- failure %>% fisher.test()

# keeping the estimates and tables and deleting the rest

# rm(list=setdiff(ls(), c("success", "enrichment_success", "failure", "enrichment_failure")))

## saving failure and success genes====

failure_genes <- unique(eg_failure$ensid.x)

success_genes <- unique(eg_success$ensid)

genes_df <- data.frame(genes = c(failure_genes, success_genes), failure_gene = c(rep(1, length(failure_genes)), rep(0, length(success_genes))), stringsAsFactors = F)

saveRDS(genes_df, "genes_failure_vs_success.rds")

# ggplot codes:====

ggplot(eg, aes(x=reorder(outcome_trait_category, outcome_trait_category,
                         function(x)-length(x)))) + geom_bar() + theme(axis.text.x = element_text(angle = 45, vjust = 0.6))

# Plot enrichment analysis====

saveRDS(success, "trial_success_table.rds")
saveRDS(failure, "trial_failure_table.rds")
saveRDS(enrichment_success, "trial_success_estimates.rds")
saveRDS(enrichment_failure, "trial_failure_estimates.rds")


# 2 x 2 tables

knitr::kable(success, caption = "Can pQTL MR-coloc associations predict of success of clinical trials?", align = 'c')

knitr::kable(failure, caption = "Can pQTL MR-coloc associations predict of success of clinical trials?", align = 'c')

# Forest plot

base_data <- tibble(mean  = c(enrichment_success$estimate, enrichment_failure$estimate), 
                    lower = c(enrichment_success$conf.int[1], enrichment_failure$conf.int[1]),
                    upper = c(enrichment_success$conf.int[2], enrichment_failure$conf.int[2]),
                    study = c("Target reached Phase 4", "Target trial terminated due to safety issues"),
                    OR = c(paste0(round(enrichment_success$estimate, 3), " (", round(enrichment_success$conf.int[1], 3), ", ", round(enrichment_success$conf.int[2], 3), ")"), paste0(round(enrichment_failure$estimate, 3), " (", round(enrichment_failure$conf.int[1], 3), ", ", round(enrichment_failure$conf.int[2], 3), ")")))

header <- tibble(study = "Target trial status", OR = "OR (95% CI)")

# empty_row <- tibble(mean = NA_real_)

df <- bind_rows(header, base_data)

styles <- fpShapesGp(
  lines = list(
    gpar(col = "darkred"),
    gpar(col = "royalblue"),
    gpar(col = "darkred")
    ),
  box = list(
    gpar(fill = "darkred"),
    gpar(fill = "royalblue"),
    gpar(fill = "darkred")
  ) 
)

# font <- "mono"
# if (grepl("Ubuntu", Sys.info()["version"])) {
#   font <- "HersheyGothicEnglish"
# }

df %>% 
  forestplot(labeltext = c(study, OR),
             xlog = TRUE,
             is.summary=c(TRUE, FALSE, FALSE),
             hrzl_lines = list("2" = gpar(lty = 2), 
                               "1" = gpar(lwd = 1, columns = 1:3, col = "#000044")),
             xlab = "Odds ratio (95% CI) of target trial status",
             txt_gp = fpTxtGp(summary = gpar(fontface = "bold")),
             shapes_gp = styles)

             