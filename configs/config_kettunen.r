# Configurations

## Authentication

json_key <- "open-targets-genetics-fc5b6cda58e5.json" # Path to json file to authenticate google cloud access

## Location of datasets

input_gc_bucket <- "gs://genetics-portal-ukbb-mr-eur/protein-gwas/KETTUNEN/harmonised" # Path to input files in a google bucket

output_gc_bucket <- "gs://genetics-portal-ukbb-mr-eur/protein-gwas/KETTUNEN/pan_filtered/" # Path to google bucket where exposure output files will be stored

outcome_gc_bucket <- "gs://genetics-portal-ukbb-mr-eur/protein-gwas/KETTUNEN/pan_outcomes_gsmr_filtered/" # store the filtered exposure-SNP matched outcomes here

## File formats

format_suffix <- TRUE # If files in input_gc_bucket have a format suffix, e.g. ‘.tsv’

parq <- FALSE # if inputs files are parquet files

## Schema

schema <- "autodetect" # Path to schema file (as an rds object) or autodetect

schema_file <- "" # if schema not specified, specify path to a sample file (representative of files in input_gc_bucket) and the pipeline will extract schema

delimiter <- "tab"

## Defining schema columns to use

schema_chr <- "hm_chrom" # Name of chromosome column in schema

schema_pos <- "hm_pos" # Name of chromosomal position column in schema (build 38)

schema_ref_allele <- "hm_other_allele" # Name of reference allele column in schema

schema_alt_allele <- "hm_effect_allele" # Name of alternative allele column in schema

schema_effect_allele_freq <- "hm_effect_allele_frequency" # Name of effect allele (typically the alt allele) frequency column in schema (Note: If eaf not provided, will use ukb eaf to annotate variant eaf)

schema_effect <- "hm_beta" # Name of effect estimate column in schema

schema_se <- "standard_error" # Name of standard error column in schema (if absent, use z-statistic and beta)

schema_p <- "p_value" # Name of p-value column in schema

schema_n <- "n_samples" # Name of sample size column in schema (for case-control studies, this will be the total sample size) (Note: one can also provide the sample size directly if they are uniform across tables)

## GC and BQ project information to use

project_id <- "open-targets-ukbb" # Google cloud project ID where input data is stored

billing_project <- "open-targets-ukbb" # Google cloud project ID to bill

dataset_name <- "kettunen" # Name of dataset that will be created in BigQuery

# filter_pval <- TRUE # create a pvalue filtered exposure and matched outcome dataset

p <- 1e-5 # specify p-value to filter on

unique_var_ds <- "mohd_hypothesis2" # any non-exposure, non-outcome stuff like list of unique variants from exposure to curate outcome datasets, goes in this dataset

outcome_ds <- "outcomes_gsmr" # name of the outcome dataset that is already gsmr formatted to match with unique_var list to curate