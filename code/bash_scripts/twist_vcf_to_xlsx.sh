#!/usr/bin/env bash

echo path to folders:
read path
echo run number:
read run

module load gatk

mkdir -p $path/tables/$run/xlsx
cd $path/annotation/$run

for vcffile in *.vcf

	do
	gunzip $vcffile

	gatk VariantsToTable -V $vcffile -O $path/tables/$run/$(echo $vcffile | cut -f 1 -d ".").table \
     	-F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO -F FORMAT -F TYPE -F RepRegion -F DP -F UMT -F VMT -F VMF \
	-F ANN -F LOF -F NMD -F CDA -F OTH -F WTD -F dbSNPBuildID -F SLO -F NSF -F R3 -F R5 -F NSN -F NSM -F G5A -F COMMON -F RS -F RV -F TPA -F CFL -F CNO -F GNO -F VLD -F ASP \
	-F ASS -F REF -F U3 -F U5 -F WGT -F MTP -F LSD -F NOC -F DSS -F SYN -F KGPhase3 -F CAF -F VC -F MUT -F KGPhase1 -F NOV -F VP -F SAO -F GENEINFO -F INT -F G5 -F OM -F PMC -F SSR \
	-F RSPOS -F HD -F PM -F CDS -F AA -F GENE -F CNT -F STRAND -F CLNSIG -F CLNALLE -F CLNORIGIN -F CLNSRC -F CLNREVSTAT -F CLNDSDB -F CLNACC -F CLNDBN -F CLNDSDBID -F CLNHGVS \
	-F CLNSRCID -F KGPilot123 -F OTHERKG -F KGValidated -F KGPROD -F PH3 -F DBVARID -F ALLELEID -F CLNVCSO -F CLNDNINCL -F ORIGIN -F MC -F CLNDN -F CLNVC -F CLNVI -F AF_EXAC \
	-F AF_ESP -F CLNSIGINCL -F CLNDISDB -F CLNDISDBINCL -F AF_TGP -F CLNSIGCONF -F CSQ \
        -GF AD -GF DP -GF GQ -GF GT -GF VF

	gzip $vcffile
done

module unload gatk

cd $path/tables/$run

for table in *.table
	do
	Rscript='/disk/work/shared/tools/R/bin/Rscript'
	t2xlsx='/disk/work/users/gb1/annotater/vcf-annotate/table_to_xlsx_QAML_run_29.R'
	$Rscript $t2xlsx
done
