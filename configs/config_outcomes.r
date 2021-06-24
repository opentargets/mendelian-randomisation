# Configurations

## Authentication

json_key <- "" # Path to json file to authenticate google cloud access

## Location of dataset

input_gc_bucket <- "" # Path to input files in a google bucket ("gs://" URL)

output_gc_bucket <- "" # Path to google bucket where exposure output files will be stored ("gs://" URL)

## File formats

format_suffix <- TRUE # If files in input_gc_bucket have a format suffix, e.g. ‘.tsv’

parq <- FALSE # if inputs files are parquet files

## Schema

schema <- "autodetect" # Path to schema file (as an rds object) or autodetect

schema_file <- "" # if schema not specified, specify path to a sample file (representative of files in input_gc_bucket) and the pipeline will extract schema

delimiter <- "tab"

## Defining schema columns to use

schema_chr <- "" # Name of chromosome column in schema

schema_pos <- "" # Name of chromosomal position column in schema (build 38)

schema_ref_allele <- "" # Name of reference allele column in schema

schema_alt_allele <- "" # Name of alternative allele column in schema

schema_effect_allele_freq <- "" # Name of effect allele (typically the alt allele) frequency column in schema (Note: If eaf not provided, will use ukb eaf to annotate variant eaf)

schema_effect <- "" # Name of effect estimate column in schema

schema_se <- "" # Name of standard error column in schema (if absent, use z-statistic and beta)

schema_p <- "" # Name of p-value column in schema

schema_n <- "" # Name of sample size column in schema (for case-control studies, this will be the total sample size) (Note: one can also provide the sample size directly if they are uniform across tables)

## GC and BQ project information to use

project_id <- "" # Google cloud project ID where input data is stored

billing_project <- "" # Google cloud project ID to bill

dataset_name <- "outcomes" # Name of dataset that will be created in BigQuery