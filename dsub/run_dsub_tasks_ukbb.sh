#!/bin/bash

# gs://genetics-portal-ukbb-mr/gcta64

# ref panel gs://open-targets-ukbb/genotypes/ukb_v3_downsampled10k
# proteins gs//genetics-portal-ukbb-mr/protein-gwas/SUN/dsub/sun/
# outcomes gs://genetics-portal-ukbb-mr/outcomes/

# image debian-10

export GOOGLE_APPLICATION_CREDENTIALS="open-targets-genetics-fc5b6cda58e5.json"

project=open-targets-ukbb
bucket_mount="gs://genetics-portal-ukbb-mr-eur"
bucket="${bucket_mount}"
logs_dir="${bucket}/pan-gsmr-outputs/logs"
# ref_dir="gs://open-targets-ukbb/genotypes/ukb_v3_downsampled10k"
# proteins_dir="${bucket}/protein-gwas/SUN/dsub/sun"
# outcomes_dir="${bucket}/outcomes_curated"
# out_dir="${bucket}/output"

gcta="gs://genetics-portal-ukbb-mr-eur/gcta64"

# trying with mount
# --input-recursive PROTEINS=$proteins_dir \

dsub \
    --provider google-cls-v2 \
    --project $project \
    --regions europe-west1 \
    --logging $logs_dir \
    --boot-disk-size 500 \
    --disk-size 3000 \
    --timeout '7d' \
    --machine-type n1-highmem-32 \
    --image gcr.io/cloud-marketplace/google/debian10 \
    --script process_gcta_batch_ukbb.sh \
    --credentials-file open-targets-genetics-fc5b6cda58e5.json \
    --tasks protein_tasks2.tsv \
    --wait