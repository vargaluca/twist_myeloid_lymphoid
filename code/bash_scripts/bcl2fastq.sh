#!/usr/bin/env bash

echo path to folders:
read path
echo run number:
read run

mkdir -p $path/raw_reads/$run/undetermined

/disk/work/shared/tools/bcl2fastq/bin/bcl2fastq \
        --input-dir $path/bcl_files/$run/Data/Intensities/BaseCalls \
        --runfolder-dir $path/bcl_files/$run \
        --output-dir $path/raw_reads/$run \
        --sample-sheet $path/bcl_files/$run/aml_all*SampleSheet_v2.csv

mv $path/raw_reads/$run/Undetermined* $path/raw_reads/$run/undetermined/

echo fastq generated for $run on $(date)
