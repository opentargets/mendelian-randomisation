

# FUNCTION

process_tables_bq <- function () {
  
  # require(bigrquery)
  
  # bq_auth()
  
  
  if(schema[[grep(schema_p, schema)]]$type=="FLOAT") { # assume beta and se are also float
    
    sapply(seq(tbl2), function (x) {
      
      tbl <- tbl2[x]
      
      print(tbl)
      
      qry <- paste0("WITH
  BIG_QUERY AS (
  SELECT
    *
  FROM (
    WITH
      BIG_SUBQUERY AS (
      SELECT
        CONCAT(", schema_chr,",':',", schema_pos, ",':',", schema_ref_allele, ",':',", schema_alt_allele,") AS SNP,",
                    schema_alt_allele," AS A1,",
                    schema_ref_allele," AS A2,
        CAST(g.string_field_4 AS numeric) AS freq,
                    SAFE_CAST(", schema_effect," AS numeric) AS beta,
                    SAFE_CAST(", schema_se," AS numeric) AS se,",
                    schema_p," AS p,
        CASE ", schema_p,"
          WHEN 0.0 THEN 1.0E-323
        ELSE ",
                    schema_p,"
      END
        AS p2,",
                    schema_n," AS N,
      FROM
        `", dataset_name,".", tbl, "` m
      INNER JOIN
        `mohd_hypothesis2.ukb_freq_dedup` g
      ON
        CONCAT(", schema_chr,",':',", schema_pos, ",':',", schema_ref_allele, ",':',", schema_alt_allele,") = g.string_field_1)
    SELECT
      SNP,
      A1,
      A2,
      freq,
      beta,
      se,
      p2 AS p,
      N,
      ROW_NUMBER() OVER (PARTITION BY SNP) row_number
    FROM
      BIG_SUBQUERY
    WHERE
      NOT (se = 0 OR beta = 0))
  WHERE
    row_number = 1)
SELECT
  SNP,
  A1,
  A2,
  freq,
  beta,
  se,
  p,
  N
FROM
  BIG_QUERY
")
      
      # bq_perform_query(qry, 
      #                  billing = billing_project, 
      #                  destination_table = bqt[x], 
      #                  write_disposition = "WRITE_TRUNCATE",
      #                  priority = "BATCH")
      
      qry <- gsub("'", '"', qry)
      
      cmd <- paste0('bq query --replace --batch --destination_table ', dataset_name2, '.', tbl, ' --use_legacy_sql=false ',  "'", qry, "'")
      
      
      system(cmd)
      
      # print("tables updated after ukb af cross-check")
      
      
    })
    
  } else {

  sapply(seq(tbl2), function (x) {
  
  tbl <- tbl2[x]
  
  print(tbl)
  
  qry <- paste0("WITH
  BIG_QUERY AS (
  SELECT
    *
  FROM (
    WITH
      BIG_SUBQUERY AS (
      SELECT
        CONCAT(", schema_chr,",':',", schema_pos, ",':',", schema_ref_allele, ",':',", schema_alt_allele,") AS SNP,",
                schema_alt_allele," AS A1,",
                schema_ref_allele," AS A2,
        CAST(g.string_field_4 AS numeric) AS freq,",
                schema_effect," AS beta,",
                schema_se," AS se,",
                schema_p," AS p,
        CASE ", schema_p,"
          WHEN '0.0' THEN '1.0E-323'
        ELSE ",
                schema_p,"
      END
        AS p2,",
                schema_n," AS N,
      FROM
        `", dataset_name,".", tbl, "` m
      INNER JOIN
        `mohd_hypothesis2.ukb_freq_dedup` g
      ON
        CONCAT(", schema_chr,",':',", schema_pos, ",':',", schema_ref_allele, ",':',", schema_alt_allele,") = g.string_field_1)
    SELECT
      SNP,
      A1,
      A2,
      freq,
      beta,
      se,
      p2 AS p,
      N,
      ROW_NUMBER() OVER (PARTITION BY SNP) row_number
    FROM
      BIG_SUBQUERY
    WHERE
      NOT (se = 'NA' OR beta = 'NA'))
  WHERE
    row_number = 1)
SELECT
  SNP,
  A1,
  A2,
  freq,
  beta,
  se,
  p,
  N
FROM
  BIG_QUERY
")
  
  # bq_perform_query(qry, 
  #                  billing = billing_project, 
  #                  destination_table = bqt[x], 
  #                  write_disposition = "WRITE_TRUNCATE", 
  #                  priority = "BATCH")
  
  qry <- gsub("'", '"', qry)
  
  cmd <- paste0('bq query --replace --batch --destination_table ', dataset_name2, '.', tbl, ' --use_legacy_sql=false ',  "'", qry, "'")
  
  system(cmd)
  
  # print("tables updated after ukb af cross-check")
  
  
})
    
  }
  
}
