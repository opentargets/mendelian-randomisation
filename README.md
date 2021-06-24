# Open Targets Genetics Mendelian randomisation pipeline

This is the github repo of the Mendelian randomisation (MR) pipeline implemented for Open Targets Genetics.

We perform two-sample pan-MR (i.e. using genome-wide significant instruments from across the genome) analysis investigating the association of genetically predicted molecular traits (proteins and metabolites) with outcomes using the [GCTA implementation of GSMR](https://cnsgenomics.com/software/gcta/#GSMR) (Generalised Summary-data-based Mendelian Randomisation).

We use most of the default settings of GSMR analysis (shown below are among the many settings, the full list of settings can be found in their [website](https://cnsgenomics.com/software/gcta/#GSMR)):

1. `--gwas-thresh` 5e-8
2. `--clump-r2` 0.05
3. `--heidi-thresh` 0.01
4. `--gsmr-snp-min` 10 (this is relaxed to 1 when datasets have limited MR results based on the min. 10 instruments req)
5. `--diff-freq` 0.5
6. `--gsmr-direction` 0 (forward MR analysis)

## Requirements

1. PLINK formatted reference genotype files split by chromosome. We use [UKBB genotype data](https://console.cloud.google.com/storage/browser/open-targets-ukbb) (not public) downsampled to 10K. UKBB genotype data were [lifted over](https://github.com/opentargets/genetics-backend/tree/master/reference_data/uk_biobank_v3) from build 37 to build 38.
2. Google Cloud Project [service account json key](https://cloud.google.com/iam/docs/creating-managing-service-account-keys)
3. Linux tools
    * [gcta](https://cnsgenomics.com/software/gcta/#Download) (`gs://genetics-portal-ukbb-mr-eur/gcta64`)
    * [bq](https://cloud.google.com/sdk/docs/install) (If using google VM, this is pre-installed as part of Google Cloud SDK)
    * [dsub](https://github.com/DataBiosphere/dsub)
    * Base R (`sudo apt-get install r-base`)
    * [RStudio](https://www.rstudio.com/products/rstudio/download-server/debian-ubuntu/)
    * Other linux-specific libraries (`sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev`)
4. R packages
    * [bigQueryR](https://code.markedmondson.me/bigQueryR/)
    * [bigrquery](https://github.com/r-dbi/bigrquery)
    * [googleCloudStorageR](https://github.com/cloudyr/googleCloudStorageR)
    * [data.table](https://github.com/Rdatatable/data.table)
    * [dplyr](https://dplyr.tidyverse.org/)
    * [arrow](https://arrow.apache.org/docs/r/) (only required if files are in parquet formar)
    
## Set up

If working on a google virtual machine (VM) instance (recommended), creating the VM with a start-up script (`gs://genetics-portal-dev-mr/start_up_script.sh`) should install all the required tools and packages:

    gcloud beta compute --project=open-targets-genetics-dev instances create gsmr \
    --zone=europe-west1-d \ 
    --machine-type=e2-medium \
    --metadata=startup-script-url=gs://genetics-portal-dev-mr/start_up_script.sh \                                     --scopes=https://www.googleapis.com/auth/cloud-platform \
    --image=debian-10-buster-v20210609 --image-project=debian-cloud \
    --boot-disk-size=100GB --boot-disk-type=pd-balanced --boot-disk-device-name=gsmr

The installations take a while and can be monitored using

    cat /var/log/daemon.log

The start-up script also installs an RStudio server. To access the RStudio server, you need to set a username and password:
    
    sudo adduser <username>
    sudo passwd <username>
    
Once these are set, exit the VM and log back in:

    gcloud compute --project "open-targets-genetics-dev" ssh --zone "europe-west1-d" "gsmr" -- -L 8787:localhost:8787
    
Enter RStudio by pointing browser to `http://localhost:8787/`

Finally,

    git clone https://github.com/opentargets/mendelian-randomisation.git

## Prepare summary statistics

Summary statistics need to be in the [cojo](https://cnsgenomics.com/software/gcta/#COJO) format.

We use harmonised summary statistics (build 38) that are harmonised using the pipelines:

- [https://github.com/EBISPOT/gwas-sumstats-harmoniser](https://github.com/EBISPOT/gwas-sumstats-harmoniser)

- [https://github.com/opentargets/genetics-sumstat-harmoniser](https://github.com/opentargets/genetics-sumstat-harmoniser)

By harmonisation, we mean that:

- the reference and alternative alleles match the corresponding alleles in the forward strand of the reference genome.
- palindromic variants, if present, are assessed using a strand consensus approach: if <u>></u>99% of the non-palindromic variants were on the forward strand, we assumed the palindromic variant would also be on the forward strand, otherwise they would be excluded from analyses.


### Outcome traits

We select outcome traits based on existing Open Targets Genetics harmonised datasets that were reported to have at least one genome-wide significant locus (source: [open-targets-genetics:200201.studies](https://console.cloud.google.com/bigquery?p=open-targets-genetics&page=table&t=studies&d=200201)(not public)). These include outcomes from Neale (round 2) and SAIGE analysis of UK Biobank traits, and outcomes from GWAS catalog. 

We format the outcome datasets first and store in BigQuery (BQ).

We use the outcome datasets in BigQuery to generate exposure-specific outcome traits (i.e. only retain genome-wide significant variants in outcome that match with exposure data) when formatting molecular traits (below). This substantially reduces the computational load and time required for GSMR analysis.

Pipeline parameters to prepare outcome sumstats are specified in the analysis config file: `configs/config_outcomes.r`

The location of the input datasets and the location where output datasets are stored must be google cloud storage (GCS) URLs.

The formatting pipeline (`prep/gsmr_format_pan_updated_outcomes.r`) uses the options set in `configs/config_outcomes.r` to ingest harmonised sumstats in BQ and uses BQ to:

- filter out rows with duplicate variant IDs
- filter out rows with null beta or standard errors
- annotate the rows with EAF from UKBB, where effect allele frequency (EAF) not provided or not available,
- filter out rows with EAF difference of more than 0.5 from UKBB EAF, where EAF of the dataset is available and used.
- annotate pval = 0 with the lowest floating point number in BQ (10<sup>-323</sup>).
- export filtered datasets to google cloud buckets (if files <u>></u> 1 GB, BQ splits the files which is then compose back and splits deleted by `prep/compose_and_delete_splits.r`)

Two files are generated in BQ:

- outcomes (harmonised sumstats)
- outcomes\_gsmr (filtered and annotated harmonised sumstats in cojo format) (this is exported to GCS)

### Molecular traits

We use the following molecular trait datasets (arranged by year):

1. ~~[LOTTA (2021)](https://www.nature.com/articles/s41588-020-00751-5)~~ (no beta/se provided)
2. [PIETZNER (2020)](https://www.nature.com/articles/s41467-020-19996-z)
3. [SCALLOP (2020)](https://www.nature.com/articles/s42255-020-00287-2)
4. [HILLARY (2019)](https://www.nature.com/articles/s42255-020-00287-2)
5. [SUN (2018)](https://www.nature.com/articles/s41586-018-0175-2)
6. [SUHRE (2017)](https://www.nature.com/articles/ncomms14357)
7. [FOLKERSEN (2017)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006706)
8. [OLLI (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5223028/)
9. [KETTUNEN (2016)](https://www.nature.com/articles/ncomms11122)
10. [DRAISMA (2015)](https://www.nature.com/articles/ncomms8208)
11. [SHIN (2014)](https://www.nature.com/articles/ng.2982)

Pipeline parameters for molecular traits are specified in the analysis config file: `configs/config_exposures.r`

The location of the input datasets and the location where output and curated outcome datasets are stored must be GCS URLs.

The formatting pipeline (`prep/gsmr_format_pan_updated_exposures.r`) uses the options set in `configs/config_exposures.r` to ingest harmonised sumstats in BQ and uses BQ to filter and annotate as for outcome traits (see above), but also:

- creates a unique p-value filtered variant list from exposure dataset (`prep/create_unique_var_list.r`)
- creates a bespoke outcome dataset for exposure based on the unique variant list (`prep/curate_outcomes.r`)

Four files are generated in BQ

- harmonised sumstats of molecular trait (&#39;mt&#39;)
- mt\_gsmr
- mt\_gsmr\_filtered
- outcomes\_mt\_gsmr\_filtered

Two of these files (mt\_gsmr\_filtered and outcomes\_mt\_gsmr\_filtered) are exported to GCS.

To run analyses, the config files are specified inside the corresponding `gsmr_format_pan_updated` R script. 

To process outcome sumstats:

    Rscript prep/gsmr_format_pan_updated_outcomes.r
    
To process exposure sumstats:

    Rscript prep/gsmr_format_pan_updated_exposure.r


## Run GSMR

Once exposure and outcome datasets are prepared, GSMR can be run in one of two ways.

1. If there are limited numbers of exposure-outcome pairs, you can produce the text files (containing location of reference, exposure, and outcome datasets), start a single google VM (if not using a local computer) and simply run the GSMR command directly on these as shown in their website. Alternatively, you can use the [GSMR R package](https://cnsgenomics.com/software/gsmr/) to perform the same analysis.

2. If there are large numbers of exposure outcome pairs and several large datasets, you can use `dsub` to use multiple google VMs and parallelise GSMR analysis.

`dsub` is an [open-source](https://github.com/DataBiosphere/dsub) command-line tool that makes it easy to submit and run batch scripts in the cloud.

We created three files that are required for `dsub` to run GSMR analysis at scale

- `dsub/tasks.tsv` - each row contains location of the exposure, outcome, reference, and output files and represents work for a single VM.
- `dsub/process_gcta_batch_ukbb.sh` - creates paths to exposure, outcome, reference, and output files and feeds to gsmr code
- `dsub/run_dsub_tasks_ukbb.sh` - `dsub` settings (e.g. VM config, disk space etc)

Once these files are prepared, we start VMs and run analyses in parallel by:

    sh dsub/run_dsub_tasks_ukbb.sh

## Extracting MR results

MR results are extracted using the scripts `extract_mr/extract_mr_results.r` and `extract_mr/extract_mr_results_per_SNP.r`.

Please note, both these scripts use a tweaked version of `extract_mr/gsmr_plot_v2.r` (original [gsmr_plot.r](https://cnsgenomics.com/software/gcta/res/gsmr\_plot.r)), which enables extraction of SNPs for an MR association with less than 10 SNPs (but only if the MR association has at least one exposure-outcome association with <u>></u> 10 instruments).

## Dataset locations (not currently publicly accessible)

- Formatted datasets: `mr_results/Location_formatted_datasets.xlsx`

- MR-analysed datasets: `mr_results/Location_MR_analysed_datasets.xls`
