#!/usr/bin/env bash

set -euo pipefail


# This is a test script to test the dry run functionality of the


test_dir=".test/wd"


rm -rf "$test_dir"


if [ ! -d "$test_dir" ]; then

mkdir -p "$test_dir"

for sample in sample1 sample2 sample3; do

    qc_folder="$test_dir/$sample/sequence_quality_control/"
    mkdir -p "$qc_folder"
    for fraction in R1 R2 se ;
    do
        touch "$qc_folder/${sample}_QC_${fraction}.fastq.gz"
    done
done
fi

set -x


snakemake -d $test_dir  --dryrun $@