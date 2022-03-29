# Using BigQuery (genomics-public-data:linkage_disequilibrium_1000G_phase_3.super_pop_EUR & open-targets-genetics:210608.variants)  for LD matrix construction

# updated version feb2022: not using v2g to assign the top gene - using vep instead

library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID 

qry1 <- paste(c("WITH bqdf4 AS (
WITH bqdf3 AS (
WITH bqdf2 AS (
WITH bqdf AS (
SELECT
    df3.*,
 CONCAT(df2.chr_id, '_', df2.position, '_', df2.ref_allele, '_', df2.alt_allele) AS varid_b38_index,
 CONCAT(df2.chr_id_b37, '_', df2.position_b37, '_', df2.ref_allele, '_', df2.alt_allele) AS varid_b37,
 df2.rs_id as rsid_index,
 df2.most_severe_consequence as most_severe_consequence_indexvar
FROM
    `open-targets-genetics.210608.variants` df2
FULL JOIN
     `open-targets-genetics-dev.linkage_disequilibrium_1000G_phase_3.super_pop_EUR` df3
ON
    CONCAT(df2.chr_id_b37, '_', df2.position_b37, '_', df2.ref_allele, '_', df2.alt_allele)=CONCAT(df3.qchrom, '_', df3.qend, '_', df3.qzeroallele, '_', df3.qoneallele)
    )
SELECT 
    bqdf.*,
    df3.most_severe_consequence as most_severe_consequence_tagvar,
    CONCAT(df3.chr_id, '_', df3.position, '_', df3.ref_allele, '_', df3.alt_allele) AS varid_b38_tag,
    rs_id AS rsid_tag
FROM
    bqdf 
LEFT JOIN
    `open-targets-genetics.210608.variants` df3
ON 
    CONCAT(df3.chr_id_b37, '_', df3.position_b37, '_', df3.ref_allele, '_', df3.alt_allele)=CONCAT(bqdf.tchrom, '_', bqdf.tend, '_', bqdf.tzeroallele, '_', bqdf.toneallele)
    )
SELECT 
    bqdf2.*,
    df4.gene_id AS gene_id_index,
    df4.fpred_max_label AS fpred_max_label_index,
FROM
    bqdf2
LEFT JOIN 
    `open-targets-ukbb.v2g.v2g_scores` df4
ON 
    CONCAT(chr_id, '_', position, '_', ref_allele, '_', alt_allele)=bqdf2.varid_b38_index
)
SELECT 
    bqdf3.*,
    df5.gene_id AS gene_id_tag,
    df5.fpred_max_label AS fpred_max_label_tag
FROM
    bqdf3
LEFT JOIN 
    `open-targets-ukbb.v2g.v2g_scores` df5
ON 
    CONCAT(chr_id, '_', position, '_', ref_allele, '_', alt_allele)=bqdf3.varid_b38_tag
WHERE 
    bqdf3.varid_b38_index IN (", toString(shQuote(snps)) ,") AND (df5.fpred_max_label IS NOT NULL OR bqdf3.fpred_max_label_index IS NOT NULL)
)
SELECT DISTINCT 
    varid_b38_index,
    rsid_index,
    fpred_max_label_index,
    gene_id_index,
    fpred_max_label_tag,
    gene_id_tag
FROM
    bqdf4
WHERE 
    POW(corr, 2) >= 0.6 OR fpred_max_label_index IS NOT NULL"), collapse = "")

snpsdf <- bq_project_query(project, qry1) %>% bq_table_download(page_size = 50000) # toggle page_size if table is too large to be parsed

# df4$r2 <- df4$corr^2

# pav <- c("missense_variant", "stop_gained", "stop_lost", "start_lost", " transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "frameshift_variant", "transcript_amplification", "inframe_insertion", "inframe_deletion", "protein_altering_variant", "TFBS_ablation") # https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html

# df4$ld_pav <- with(df4, ifelse(most_severe_consequence %in% pav & r2 >= 0.6, "Yes", "No")) # NEED TO FIX THIS

# df4 <- df4[which(!duplicated(df4$varid_b38_index)),]

# snpsdf <- df4 %>% select(varid_b38, rs_id, most_severe_consequence, gene_id_topv2g, ld_pav)
