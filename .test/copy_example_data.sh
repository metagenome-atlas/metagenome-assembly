#!/usr/bin/env bash

set -euo pipefail

git clone https://github.com/metagenome-atlas/example_data.git

src_dir="example_data/reads/stub"

test_dir="$1"

mkdir -p "$test_dir"/..
mkdir "$test_dir"



#copy files
for sample in Mycoplasma Streptococcus; do

    qc_folder="$test_dir/$sample/sequence_quality_control/"
    mkdir -p "$qc_folder"
    for fraction in R1 R2 ;
    do
        mv $src_dir/${sample}_${fraction}.fastq.gz "$qc_folder/${sample}_QC_${fraction}.fastq.gz"
    done
done
