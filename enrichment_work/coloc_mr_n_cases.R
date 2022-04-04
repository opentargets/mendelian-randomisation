# Value of pan-MR beyond cis-coloc

library(dplyr)
library(ggplot2)

df <- readRDS("~/mr_prot_unfiltered_dataset_v1_v2_without_egger.rds")

df_sig <- df %>% filter(bxy_pval < 0.0005 | coloc_h4 > 0.8)

df_sig_gsmr_coloc <- df_sig %>% filter(!is.na(exp_out_gsmr))

eg <- df_sig_gsmr_coloc %>% tidyr::unchop(cols = c("outcome_trait_efo"))

eg$key <- with(eg, paste0(ensid, "_", outcome_trait_efo))

eg2 <- eg %>% group_by(key) %>% filter(!is.na(exp_out_coloc) & coloc_h4 > 0.8) %>% select(key) %>% distinct()

lst <- lapply(seq(eg2$key), function (x) {
  
  eg3 <- eg[which(eg$key==eg2$key[x]),]
  
  cond <- anyNA(eg3$exp_out_coloc)
  
  if (cond==TRUE) {
    
    lst <- eg3
    
  }
  
})

lst2 <- do.call(rbind, lst)

ukeys <- unique(lst2$key) # 370 target-trait associations from different studies had MR association but only some of the studies had coloc

eg3 <- lst2[which(lst2$key==ukeys[5]),]

eg3 %>% ggplot(aes(x=reorder(outcome, n_cases, FUN = sum), y=n_cases, fill = factor(coloc))) +
  geom_bar(stat='identity') + theme(axis.text.x=element_text(angle = 45, hjust = 1)) + labs(title=paste0(unique(eg3$hgnc_protein),  " - ", unique(eg3$outcome_trait_efo)), x="study_id") + geom_text(aes(label=n_cases), position=position_dodge(width=0.9), vjust=-0.25)
                                    
                                                                                            


# lst2$key2 <- with(lst2, paste0(key, "_", coloc))
# 
# lst3 <- lst2 %>%
#   group_by(key2) %>%
#   summarise(mean = mean(n_cases), n = n()) %>% filter(!is.na(mean))
# 
# coloc <- lst3[grep("Yes", lst3$key2),]
# 
# no_coloc <- lst3[grep("No", lst3$key2),]


# rm(list=setdiff(ls(), "lst2"))

# other way (none)

# df <- readRDS("~/mr_prot_unfiltered_dataset_v1_v2_without_egger.rds")
# 
# df_sig <- df %>% filter(bxy_pval < 0.0005 | coloc_h4 > 0.8)
# 
# df_sig_gsmr_coloc <- df_sig %>% filter(!is.na(exp_out_gsmr))
# 
# eg <- df_sig_gsmr_coloc %>% tidyr::unchop(cols = c("outcome_trait_efo"))
# 
# eg$key <- with(eg, paste0(ensid, "_", outcome_trait_efo))
# 
# eg2 <- eg %>% group_by(key) %>% filter(!is.na(exp_out_gsmr)) %>% select(key) %>% distinct()
# 
# lst <- lapply(seq(eg2$key), function (x) {
#   
#   eg3 <- eg[which(eg$key==eg2$key[x]),]
#   
#   cond <- anyNA(eg3$exp_out_gsmr)
#   
#   if (cond==TRUE) {
#     
#     lst <- eg3
#     
#   }
#   
# })
# 
# lst2 <- do.call(rbind, lst) # none

