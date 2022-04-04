# Compare the overall direct score of all target-trait pairs with MR evidence with the overall direct score of all target-traits without any genetic evidence

library(dplyr)

# Load target-trait pairs supported by MR

df_mr_filtered <- data.table::fread("~/prot_mr_filtered.tsv", stringsAsFactors = F)

df_mr_filtered$outcome_trait_efo[df_mr_filtered$outcome_trait_efo==""] <- NA

df_mr <- df_mr_filtered %>% tidyr::separate_rows(outcome_trait_efo, sep = ",") %>% filter(!is.na(outcome_trait_efo) & !is.na(exp_out_gsmr)) %>% select(ensid, outcome_trait_efo) %>% distinct()

df_mr$key <- with(df_mr, paste0(ensid, "_", outcome_trait_efo))

# reproduce Daniel's analysis====

library(bigrquery)

project <- "open-targets-genetics" # replace this with your project ID

qry <- "SELECT DISTINCT
  genetic.diseaseId as disease,
  genetic.targetId as target,
  genetic.score as score
FROM
  `open-targets-prod.platform.associationByDatasourceDirect` genetic
WHERE
  datasourceId='ot_genetics_portal'"

tbl <- bq_project_query(project, qry) %>% bq_table_download()

tbl$key <- with(tbl, paste0(target, "_", disease))

mean_tbl <- mean(tbl$score) # 0.224

# merge tbl with MR table

dftbl <- merge(tbl, df_mr, by = "key")

mean_dftbl <- mean(dftbl$score) # 0.307

# are the two means significantly different?

sigmean <- t.test(dftbl$score, tbl$score)

sigmean$p.value # 1.280928e-92, so yes

# Load target-trait pairs supported by only non-OTG evidence ==== 

library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID

qry <- "# Creating a table of target-trait pairs that don't have genetic support at all

SELECT DISTINCT
  nongenetic.diseaseId as disease,
  nongenetic.targetId as target,
  overall.score as score,
  nongenetic.datatypeId
FROM
  `open-targets-prod.platform.associationByDatasourceDirect` nongenetic
LEFT JOIN 
    `open-targets-prod.platform.associationByOverallDirect` overall
ON
    nongenetic.targetId = overall.targetId AND nongenetic.diseaseId = overall.diseaseId
WHERE
  nongenetic.datatypeId!='genetic_association'"

tbl <- bq_project_query(project, qry) %>% bq_table_download()

tbl$key <- with(tbl, paste0(target, "_", disease))

tbl2 <- tbl %>% group_by(key) %>% filter(score == max(score)) # retaining unique target-trait pairs

# counts

# tbl2 <- tbl %>% filter(!duplicated(key))

# tbl <- tbl %>% group_by(key) %>% mutate(otg = ifelse(any(datasourceId=="ot_genetics_portal"), 1, 0))

# tbl2 <- tbl %>% filter(otg==1) # 322,360 target-trait pairs with OTG evidence

#tbl2 <- tbl %>% filter(otg==0) # most don't have target-trait evidence

mean_tbl <- mean(tbl$score) # 0.05

# merge non-genetic tbl with MR table

dftbl <- merge(tbl, df_mr, by = "key") # 2036 unique target-trait pairs out of 186,4814, so only covers 0.1% of target-trair pairs with non-genetical support

mean_dftbl <- mean(dftbl$score) # 0.1

# are the two means significantly different?

sigmean <- t.test(dftbl$score, tbl$score)

sigmean$p.value # 2.046804e-38, so yes


# Does MR help with target prioritization when L2G < 0.5?====

library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID

qry <- "SELECT DISTINCT
  genetic.diseaseId as disease,
  genetic.targetId as target,
  genetic.score as score
FROM
  `open-targets-prod.platform.associationByDatasourceDirect` genetic
WHERE
  datasourceId='ot_genetics_portal' AND score < 0.5"

tbl <- bq_project_query(project, qry) %>% bq_table_download()

tbl$key <- with(tbl, paste0(target, "_", disease))

mean_tbl <- mean(tbl$score) # 0.173

dftbl <- merge(tbl, df_mr, by = "key")

mean_dftbl <- mean(dftbl$score) # 0.204

# are the two means significantly different?

sigmean <- t.test(dftbl$score, tbl$score)

sigmean$p.value # 1.326432e-26, so yes

nrow(dftbl)/nrow(df_mr)*100 # 10% of target-trait pairs with MR evidence overlap with OTG target-trait pairs where L2G < 0.5

# Summary

# 2036 unique target-trait pairs (0.1%) that are not supported by 'genetic_association' (but supported by other sources of evidence in the platform) are supported by MR

# The mean overall score (not L2G score) of the 2036 target-trait pairs is 0.1. When compared to the mean overall score of the all non-genetically supported target-trait pairs - 0.05, it is significant. 

# Given the value of pQTL MR in causal gene prioritization, pQTL MR can be help guide effector gene identification where L2G < 0.5 (i.e. for a given target-trait pair linked via a GWAS locus, the probability of the target being part of the gold standard set is less than 50%).

# 270,968 out of 306,660 target-trait pairs have L2G < 0.5

# There are 2542 target-trait pairs (10% of MR target-trait pairs) out of 270,968 target-trait pairs (~1% of L2G < 0.5 target-trait pairs) where L2G < 0.5 and where pQTL MR can help prioritise the right effector gene.

# The mean L2G score of the overlapping dataset with MR evidence is 0.204. When compared to the mean L2G score of the L2G < 0.5 dataset - 0.173, it is significant



# Conclusion

- The above analyses suggest that there is some value of MR (albeit small given the low proportions of target-trait pairs that overlap - mainly limited by what targets are assayed in the protein GWA studies) beyond what is simply supported by genetic colocalisation as evidenced by the relatively higher L2G and overall scores for target-trait supported by MR irrespective of genetic colocalisation.  

- When L2G score for a given target-trait pair is high (e.g. > 0.5), pQTL MR has limited value. However, pQTL MR may be of particularly use to resolve the likely causal gene for traits with L2G < 0.5.  

- pQTL MR may also help provide genetic support to additional target-trait pairs that are presently not genetically supported in the platform, but this remains a very small proportion. 