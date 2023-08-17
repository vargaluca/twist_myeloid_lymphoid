#!/usr/bin/env bash

echo path to folders:
read path
echo run number:
read run

mkdir -p $path/tables/$run
cd $path/maf_files/$run

Rscript='/disk/work/shared/tools/R/bin/Rscript'
xlsx='/disk/work/users/lv1/twist_myeloid_lymphoid/code/R'

for maffile in *.maf.gz
	do
	name=$(echo $maffile | cut -f 1 -d ".")

#	zcat $maffile > $path/tables/$run/$(echo $maffile | cut -f 1 -d ".").table
	$Rscript $xlsx/table_to_xlsx.R $path/maf_files/$run/$maffile $path/tables/$run/$name.xlsx

	if [ -f $path/tables/$run/$name.xlsx ]
	then
		if [ -s $path/tables/$run/$name.xlsx ]
		then
			echo "$name.xlsx successfully created"
		else
			echo "$name.xlsx exists but is empty"
		fi
	else
		echo "$name.xlsx doesn't exist"
	fi
done
