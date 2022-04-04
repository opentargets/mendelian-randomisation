# DOSE-RESPONSE PLOT

# x-axis: exposure
# y-axis: outcome

# Example: ST2_NEALE2_30150 (SCALLOP_2020: ST2 and Eosinophil count: mixed instruments)

files <- dir("~/mr_results_snps/3_annotation_eggerwm_outcometraits_ensid", full.names = T)

prot <- c("SUN_2018", "SUHRE_2017", "SCALLOP_2020", "PIETZNER_2020", "OLLI_2017", "HILLARY_2019", "FOLKERSEN_2017")

files2 <- grep(paste(prot, collapse = "|"), files, value = T)

lst2 <- lapply(seq(files2), function (x) {
  
  df <- readRDS(files2[x])
  
})

df2 <- bind_rows(lst2)

df2$exp_out <- with(df2, paste0(exposure, "_", outcome))

df3 <- df2 %>% select(exp_out, snp, bzx, bzx_se, bzx_pval, bzx_pval_mantissa, bzx_pval_exponent, bzy, bzy_se, bzy_pval, bzy_pval_mantissa, bzy_pval_exponent, ensid, cis_trans)

# Isolate the data

eg <- df3[df3$exp_out=="ST2_NEALE2_30150",]

snps <- gsub(":", "_", unique(eg$snp))

source("~/SNP_v2g_vep_annotation.R")

# Integrating the annotated set with the SNP set

df4$varid <- gsub("_", ":", df4$varid)

df5 <- merge(eg, df4, by.x = "snp", by.y = "varid", all.x = TRUE)

# Plot (scatter plot)

library(MendelianRandomization)

mr_plot(mr_input(bx = df5$bzx, bxse = df5$bzx_se, by = df5$bzy, byse = df5$bzy_se, snps = df5$snp, exposure = "IL1RL1 protein (SCALLOP_2020)", outcome = "Eosinophil count (NEALE2_30150)"), line="ivw", labels = T)
        
        
# plot (regional association plots)

library(gassocplot2)

markers <- gassocplot2::test_stack_assoc_plot_markers

z <- gassocplot2::test_stack_assoc_plot_associations

corr <- gassocplot2::test_corr # this is correlation not correlation squared and has to be ordered in the same way as the markers data frame

plot <- stack_assoc_plot(markers, z, corr, traits=c("Trait 1", "Trait 2"))

stack_assoc_plot_save(plot, "stack_assoc_plot_test.png", 2)

