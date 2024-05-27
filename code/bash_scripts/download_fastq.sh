#!/usr/bin/env bash

echo download path:
read path
echo run number:
read run
echo download project number:
read project

mkdir -p $path/raw_reads/$run/undetermined

/disk/work/shared/tools/bs/bs project download --id $project --extension fastq.gz -o $path/raw_reads/$run
mv $path/raw_reads/$run/Undetermined* $path/raw_reads/$run/undetermined/

cd $path/raw_reads/$run

	for folder in *
		do
		mv $folder/* $path/raw_reads/$run
	done

rmdir *
echo download finished on $(date)


