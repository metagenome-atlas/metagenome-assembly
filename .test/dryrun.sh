#!/usr/bin/env bash

set -euo pipefail


# This is a test script to test the dry run functionality of the


test_dir=".test/wd"


rm -rf "$test_dir"


mkdir -p "$test_dir"

for sample in sample1 sample2 sample3; do

    qc_folder="$test_dir/$sample/sequence_quality_control/"
    mkdir -p "$qc_folder"
    for fraction in R1 R2 se nanopore ;
    do
        touch "$qc_folder/${sample}_QC_${fraction}.fastq.gz"
    done
done


set -x

echo "Dryrun with metaspades"

snakemake -d $test_dir  --dryrun $@





echo "Dryrun with megahit"

snakemake -d $test_dir  --dryrun --config assembler="megahit" $@

echo "Dryrun with filtering"

snakemake -d $test_dir  --dryrun --config assembler="megahit" filter_contigs=True $@


echo "Dryrun with different options"

snakemake -d $test_dir  --dryrun --config merge_pairs_before_assembly=False $@


echo "Dryrun with nanopore reads"

snakemake -d $test_dir  --dryrun --config longread_type="nanopore" $@






# create folder for single ends
set +x
fraction="se"

rm -rf "$test_dir"
mkdir -p "$test_dir"

for sample in sample1 sample2 sample3; do

    qc_folder="$test_dir/$sample/sequence_quality_control/"
    mkdir -p "$qc_folder"

    

    touch "$qc_folder/${sample}_QC_${fraction}.fastq.gz"

done

set -x


echo "Dryrun with single end reads and megahit"

snakemake -d $test_dir  --dryrun --config assembler="megahit" paired_end=False $@



echo -e '\U0001F64C \U0001F389'