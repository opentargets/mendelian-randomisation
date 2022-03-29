library(bigrquery)
library(dplyr)

project <- "open-targets-genetics" # replace this with your project ID

qry1 <- "WITH bqdf2 AS (
WITH bqdf AS (
SELECT
  CONCAT(left_chrom, ':', left_pos, ':', left_ref, ':', left_alt) AS varid_left,
  left_study,
  right_study,
  right_phenotype,
  right_gene_id,
  coloc_n_vars,
  coloc_h1,
  coloc_h2,
  coloc_h3,
  coloc_h4,
  coloc_h4_h3,
  left_var_right_study_beta,
  left_var_right_study_se,
  left_var_right_study_pval
FROM
  `open-targets-genetics-dev.genetics_dev.variant_disease_coloc`
WHERE 
    right_type = 'pqtl'
)
SELECT 
    bqdf.*,
    df.beta as gwas_beta,
    (df.beta_ci_upper - df.beta_ci_lower)/3.92 as gwas_se, #https://handbook-5-1.cochrane.org/chapter_7/7_7_7_2_obtaining_standard_errors_from_confidence_intervals_and.htm
    log(df.odds_ratio) as gwas_beta_from_or,
    (log(oddsr_ci_upper) - log(oddsr_ci_lower))/3.92 as gwas_se_from_orci,
    pval
FROM 
    bqdf
INNER JOIN
  `open-targets-genetics-dev.genetics_dev.variant_disease` df
ON #WHY DOES JEREMY USE TAG VARIANTS FROM V2D TO JOIN WITH COLOC DATA? 
CONCAT(varid_left, ':', left_study)=CONCAT(lead_chrom, ':', lead_pos, ':', lead_ref, ':', lead_alt, ':', study_id)
)
SELECT
  *
FROM (
  SELECT
    *,
    ROW_NUMBER() OVER (PARTITION BY CONCAT(varid_left, ':', left_study, ':', right_study, ':', right_gene_id)) row_number
  FROM
    bqdf2)
WHERE
  row_number = 1"

coloc <- bq_project_query(project, qry1) %>% bq_table_download(page_size = 30000) # toggle page_size if table is too large to be parsed

coloc$row_number <- NULL