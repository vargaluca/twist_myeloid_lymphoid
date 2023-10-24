#!/usr/bin/env bash

echo download path:
read path
echo run number:
read run
echo download run id:
read run_id

mkdir -p $path/bcl_files/$run

/disk/work/shared/tools/bs/bs run download --id $run_id -o $path/bcl_files/$run

echo bcl2 download finished for $run_id on $(date)


