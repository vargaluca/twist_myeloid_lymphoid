#!/usr/bin/env bash

echo path to folders:
read path
echo run number:
read run

reference='/disk/work/shared/data/genomes/ucsc.hg19_bwa_index/ucsc.hg19.fasta.gz'
germline='/disk/work/diagnostics/DLBC/lv1/DLBCL_genclass_K21/gnomad_somatic/af-only-gnomad.raw.sites.lifted_hg19.vcf'
somatic='/disk/work/diagnostics/DLBC/lv1/DLBCL_genclass_K21/small_exac_common/small_exac_common_3.lifted_hg19.vcf'
pon='/disk/work/diagnostics/DLBC/bb1/DLBCL_genclass_K21/variants/pon/gatk_wes_panel/Mutect2-exome-panel.hg19.vcf'


## VARIANT CALLING USING GATK


mkdir -p $path/variant_calls/$run/f1r2
mkdir -p $path/variant_calls/$run/unfiltered_vcf
cd $path/recal_aligned_reads/$run

module load gatk

##running mutect2 for variant calls

	for bamfile in *.bam
                do gatk Mutect2 -R $reference -I $bamfile --germline-resource $germline --panel-of-normals $pon \
                --f1r2-tar-gz $path/variant_calls/$run/f1r2/$(echo $bamfile | cut -d "." -f 1)-f1r2.tar.gz \
                -O $path/variant_calls/$run/unfiltered_vcf/$(echo $bamfile | cut -d "." -f 1).unfiltered.vcf.gz
		echo $(echo $bamfile | cut -f 1 -d "_") mutect2 finished >> $path/log_files/$run.log.txt
                done

echo UNFILTERED VCFs FINISHED ON $(date) >> $path/log_files/$run.log.txt

cd $path/variant_calls/$run/f1r2

##filtering out artifacts from FFPE samples

mkdir -p $path/variant_calls/$run/read_orientation_model

	for filename in *tar.gz
		do gatk LearnReadOrientationModel -I $filename \
		-O $path/variant_calls/$run/read_orientation_model/$(echo $filename | sed 's/f1r2/read-orientation-model/')
		done

cd $path/recal_aligned_reads/$run
mkdir -p $path/variant_calls/$run/pileup_summaries

	for bamfile in *.bam
		do gatk GetPileupSummaries -I $bamfile -V $somatic -L $somatic \
		-O $path/variant_calls/$run/pileup_summaries/$(echo $bamfile | cut -d "." -f 1)-pileupsummaries.table
		echo $(echo $bamfile | cut -f 1 -d "_") pileupsummaries finished >> $path/log_files/$run.log.txt
		done

cd $path/variant_calls/$run/pileup_summaries
mkdir -p $path/variant_calls/$run/contamination

	for pileup in *.table
		do gatk CalculateContamination -I $pileup \
		-O $path/variant_calls/$run/contamination/$(echo $pileup | sed 's/pileupsummaries/contamination/')
		echo $(echo $bamfile | cut -f 1 -d "_") contamination finished >> $path/log_files/$run.log.txt
		done

##filtering vcf files

cd $path/variant_calls/$run/unfiltered_vcf
mkdir -p $path/variant_calls/$run/filtered_vcf

	for vcffile in *unfiltered.vcf.gz
		do gatk FilterMutectCalls -R $reference -V $vcffile \
		--contamination-table $path/variant_calls/$run/contamination/$(echo $vcffile | cut -d "." -f 1)-contamination.table \
		--ob-priors $path/variant_calls/$run/read_orientation_model/$(echo $vcffile | cut -d "." -f 1)-read-orientation-model.tar.gz \
		-O $path/variant_calls/$run/filtered_vcf/$(echo $vcffile | sed 's/unfiltered/filtered/')
		echo $(echo $bamfile | cut -f 1 -d "_") filtered vcf finished >> $path/log_files/$run.log.txt
		done

module unload gatk

echo GATK VARIANT CALLING FINISHED ON $(date)


## ANNOTATION USING VEP


vep_annot='/disk/work/shared/data/vep'

mkdir -p $path/annotation/$run
cd $path/variant_calls/$run/filtered_vcf

module load vep

	for vcffile in *vcf.gz
        	do
		outfile=$(echo $vcffile | sed 's/filtered/annotated/')

        	vep --dir_cache $vep_annot --cache --fasta $reference --offline --everything -i $vcffile \
        	--vcf --pick --force_overwrite --assembly GRCh37 --fork 8 \
        	-o $path/annotation/$run/$outfile

        	mv $path/annotation/$run/$outfile \
        	$path/annotation/$run/$(echo $outfile | sed 's/vcf.gz/vcf/')
		echo $(echo $bamfile | cut -f 1 -d "_") annotation finished >> $path/log_files/$run.log.txt
	done

module unload vep

echo VEP ANNOTATION FINISHED ON $(date)


## CREATING MAF FILES FROM ANNOTATED VCF


module load samtools
module load vcf2maf

vcf2maf='/disk/work/shared/tools/vcf2maf/vcf2maf.pl'
tabix='/work/shared/tools/ensembl-vep-release-106/htslib/tabix'
reference='/disk/work/shared/data/genomes/ucsc.hg19_bwa_index/ucsc.hg19.fasta.gz'
vep_exec='/disk/work/shared/tools/ensembl-vep'
vep_data='/disk/work/shared/data/vep'

cd $path/annotation/$run
mkdir -p $path/maf_files/$run

	for vcf in *vcf
        	do
		perl $vcf2maf --input-vcf $vcf \
        	--output-maf $path/maf_files/$run/$(echo $vcf | sed 's/annotated.vcf/maf/') \
        	--retain-info --retain-fmt \
        	--tumor-id $(echo $vcf | cut -f 1 -d "_") --inhibit-vep --tabix-exec $tabix \
        	--ref-fasta $reference --vcf-tumor-id $(echo $vcf | cut -f 1 -d "_")

        	gzip $path/maf_files/$run/$(echo $vcf | sed 's/annotated.vcf/maf/')
		echo $(echo $bamfile | cut -f 1 -d "_") maf finished >> $path/log_files/$run.log.txt
        done

module unload samtools
module unload vcf2maf

mkdir -p $path/tables/$run
cd $path/maf_files/$run

## CREATING VARIANT TABLES

Rscript='/disk/work/shared/tools/R/bin/Rscript'
xlsx='/disk/work/users/lv1/twist_myeloid_lymphoid/code/R/table_to_xlsx.R'
combined='/disk/work/users/lv1/twist_myeloid_lymphoid/code/combined_table.R'

for maffile in *.maf.gz
        do
        name=$(echo $maffile | cut -f 1 -d ".")

#       zcat $maffile > $path/tables/$run/$(echo $maffile | cut -f 1 -d ".").table
        $Rscript $xlsx $path/maf_files/$run/$maffile $path/tables/$run/$name.xlsx

        if [ -f $path/tables/$run/$name.xlsx ]
        then
                if [ -s $path/tables/$run/$name.xlsx ]
                then
                        echo "$name.xlsx successfully created" >> $path/log_files/$run.log.txt
                else
                        echo "$name.xlsx exists but is empty" >> $path/log_files/$run.log.txt
                fi
        else
                echo "$name.xlsx doesn't exist" >> $path/log_files/$run.log.txt
        fi
done

$Rscript $combined $path/maf_files/$run/ $path/tables/$run/$run-combined_table.xlsx

echo VARIANT TABLES FINISHED FOR $run ON $(date)

