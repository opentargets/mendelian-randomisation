# Preparing MR datasets to share with partners (Feb 2022 release)

library(dplyr)

dfmr2 <- readRDS("~/mr_prot_filtered_dataset_v1_v2.rds") # filtered and annotated MR dataset

dfsnp2 <- readRDS("~/mr_prot_snp_filtered_dataset_v1_v2.rds") # filtered and annotated SNP dataset

dfmr2 <- dfmr2 %>% select(exp_out_gsmr_coloc, Data, exposure, protein_trait, hgnc_protein, ensid, outcome, outcome_trait, outcome_trait_efo, nsnp, n_cases, n_initial, cis_trans_mr, coloc, varid_left, coloc_h4, coloc_h4_h3,  bxy, bxy_se, bxy_pval, bxy_pval_mantissa, bxy_pval_exponent, mr_egger, mr_egger_se, mr_egger_p, mr_egger_p_mantissa, mr_egger_p_exponent, mr_wm, mr_wm_se, mr_wm_p, mr_wm_p_mantissa, mr_wm_p_exponent, pav_cismr, fpred_max_label_index_mr, fpred_max_label_tag_mr, v2g_mr, hgnc_v2g_mr)

dfsnp2 <- dfsnp2 %>% select(exp_out_gsmr_coloc, snp, cis_trans, bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, fpred_max_label_index, fpred_max_label_tag, pav_cis, hgnc_v2g)

# convert list columns to character (because we are saving as csv file)

list_to_char1 <- which(sapply(dfmr2, class)=="list")

sapply(seq(list_to_char1), function (x) { dfmr2[,list_to_char1[x]] <<- as.character(dfmr2[,list_to_char1[x]])})

# renaming varid_left

names(dfmr2)[grep("varid_left", names(dfmr2))] <- "coloc_snp"

# save as csv.gz

write.table(dfmr2, file=gzfile("~/protein_MR_March2022.tsv.gz"), sep = "\t", quote = F, row.names = F)
write.table(dfsnp2, file=gzfile("~/protein_MR_snp_March2022.tsv.gz"), sep = "\t", quote = F, row.names = F)

# exporting to google cloud

version <- "v2"

system(paste0("gsutil cp protein_MR_March2022.tsv.gz gs://genetics-portal-dev-mr/mr_results_", version, "/partners/"))

system(paste0("gsutil cp protein_MR_snp_March2022.tsv.gz gs://genetics-portal-dev-mr/mr_results_", version, "/partners/"))

system(paste0("gsutil cp README_protein_MR_March2022.txt gs://genetics-portal-dev-mr/mr_results_", version, "/partners/"))

system(paste0("gsutil cp preparing_mr_datasets_partners.R gs://genetics-portal-dev-mr/mr_results_", version, "/partners/"))

# README

# "README_protein_MR_March2022.txt" done

# Export to otcoregen farm (instructions)

# Please provide datasets to be shared with Open Targets partners as described below.
# Data not shared in this manner may remain unavailable to the partners until corrected.
# Include 1) short, one-sentence description of the dataset for the intranet page & 2) a brief README/description of the files being shared and how they relate to one another.
# 
# Copy individual releases to sub-directories of your project.
# /lustre/scratch118/opentargets/opentargets/OTAR***/share/
#   E.g., Release #1:
# /lustre/scratch118/opentargets/opentargets/OTAR***/share/Release_1
# Place subsequent releases in different directories.
# /lustre/scratch118/opentargets/opentargets/OTAR***/share/Release_2
# /lustre/scratch118/opentargets/opentargets/OTAR***/share/Release_3
# Ensure that all the copied data is group read/writeable.
# chmod –R 770 /lustre/scratch118/opentargets/opentargets/OTAR***/share/Release_1
# 
# 
# NOTES:
#   Data not copied to a directory under /share/ may go unnoticed or completely deleted without warning, thus remain unavailable for sharing.
# 
# Data presented for sharing may be archived off the system once it is available to the partners. If you would like to retain access to the uploaded data during the project, it should be stored elsewhere.
# 
# Any questions, please speak to either Sarah Young (sy3@sanger.ac.uk) or Stuart Horswell (sh40@sanger.ac.uk) for more details.

# how I actually do it (in sanger laptop, not google VM)

# ssh farm5-login
# cd /lustre/scratch119/opentargets/opentargets/OTAR034/share