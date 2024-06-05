#!/usr/bin/env bash

echo path to folders:
read path
echo run number:
read run


## RUNNING FASTQC ON RAW READS

mkdir -p $path/tables/$run
mkdir -p $path/log_files
touch $path/log_files/$run.log.txt

mkdir -p $path/raw_fastqc/$run
cd $path/raw_reads/$run

	for filename in *.fastq.gz
	       	do
		 /disk/work/shared/tools/FastQC/fastqc -o $path/raw_fastqc/$run $filename
		echo $(echo $filename | cut -f 1 -d "_") fastqc finished >> $path/log_files/$run.log.txt
        done

cd $path/raw_fastqc/$run
multiqc . --filename $run-raw_multiqc

mv $run-raw_multiqc.html $path/tables/$run/

echo FASTQC AND MULTIQC FINISHED ON $(date) >> $path/log_files/$run.log.txt


## TRIMMING RAW READS WITH AGENT


module load agent

mkdir -p $path/trimmed_reads/$run
cd $path/raw_reads/$run

	for read1 in *R1*.fastq.gz
        	do
		agent trim -v2 -out $path/trimmed_reads/$run/$(echo $read1 | sed 's/_L001_R1_001.fastq.gz/.trimmed/') \
        	-fq1 $read1 \
        	-fq2 $(echo $read1 | sed 's/R1/R2/')
		echo $(echo $read1 | cut -f 1 -d "_") trimming finished >> $path/log_files/$run.log.txt
        done

module unload agent
echo TRIMMING FINISHED ON $(date) >> $path/log_files/$run.log.txt


## RUNNING FASTQC ON TRIMMED READS


mkdir -p $path/trimmed_fastqc/$run
cd $path/trimmed_reads/$run

	for filename in *.fastq.gz
        	do
		/disk/work/shared/tools/FastQC/fastqc -o $path/trimmed_fastqc/$run $filename
		echo $(echo $filename | cut -f 1 -d "_") trimmed fastqc finished >> $path/log_files/$run.log.txt
        done

cd $path/trimmed_fastqc/$run
multiqc . --filename $run-trimmed_multiqc
mv $run-trimmed_multiqc.html $path/tables/$run/

echo TRIMMED FASTQC FINISHED ON $(date) >> $path/log_files/$run.log.txt


## ALIGNING TRIMMED FASTQ TO REFERENCE GENOME WITH BWA-MEM


reference='/disk/work/shared/data/genomes/ucsc.hg19_bwa_index/hg19.fa.gz'
mkdir -p $path/aligned_reads/$run
cd $path/trimmed_reads/$run

	for read1 in *trimmed_R1.fastq.gz
		do
		header=$(zcat $read1 | head -n 1); id=$(echo $header | cut -f 3-4 -d ":" | sed 's/@//' | sed 's/:/./g')
        	sample=$(echo $read1 | cut -f 1 -d "_")
        	platform_unit=$(echo $header | cut -f 1,3,4 -d ":" | sed 's/@//' | sed 's/:/./g')

        	/disk/work/shared/tools/bwa/bwa mem -R $(echo "@RG\tID:$id\tPL:illumina\tPU:$platform_unit\tSM:$sample\tLB:$sample") -C -t 8 \
		$reference $read1 \
        	$(echo $read1 | sed 's/R1/R2/') | \
	       	/disk/work/shared/tools/samtools/samtools view -b > $path/aligned_reads/$run/$(echo $read1 | sed 's/trimmed_R1.fastq.gz/aligned.bam/')

		echo $(echo $read1 | cut -f 1 -d "_") alignment to hg19 finished >> $path/log_files/$run.log.txt
	done

## sorting .bam files using samtools

cd $path/aligned_reads/$run

module load samtools

	for bamfile in *.bam
        	do
        	samtools sort $bamfile -o $(echo $bamfile | cut -f 1 -d ".").sorted.bam
		echo $(echo $bamfile | cut -f 1 -d "_") sorting finished >> $path/log_files/$run.log.txt
	done

# line 100
echo BWA-MEM ALIGNMENT FINISHED ON $(date) >> $path/log_files/$run.log.txt

# creating bamqc and mutiqc for the bam files

mkdir -p $path/bam_qc/$run
cd $path/aligned_reads/$run

module load qualimap

covered='/disk/work/diagnostics/TWLM/lv1/twist_myeloid_lymphoid/bed_files/Target_bases_covered_by_probes_Amplikon_LymphoidMyeloid_1X_TE-95491860_hg19.bed'

	for bamfile in *sorted.bam
		do
        	outdir=$(echo $bamfile | cut -f 1 -d ".")
        	mkdir -p $path/bam_qc/$run/$outdir

        	qualimap bamqc -bam $bamfile -outdir $path/bam_qc/$run/$outdir -gff $covered --java-mem-size=4G

        	mv $path/bam_qc/$run/$outdir/qualimapReport.html \
        	$path/bam_qc/$run/$outdir/$(echo $bamfile | cut -f 1 -d ".").qualimapReport.html
        	echo $path/bam_qc/$run/$outdir >> $path/bam_qc/$run/file_paths.txt

		echo $(echo $bamfile | cut -f 1 -d "_") bamqc finished >> $path/log_files/$run.log.txt
	done

cd $path/bam_qc/$run

multiqc -n $path/bam_qc/$run/reports/$run.aml_all-bammultiqc --file-list file_paths.txt
mv $path/bam_qc/$run/reports/$run.aml_all-bammultiqc.html $path/tables/$run/

rm file_paths.txt
module unload qualimap
module unload samtools

echo BAMQC FINISHED ON $(date) >> $path/log_files/$run.log.txt


## DEDUPLEXING ALIGNED BAM FILES


mkdir -p $path/dedup_aligned_reads/$run
mkdir -p $path/dedup_aligned_reads/tmp/$run
cd $path/aligned_reads/$run

module load gatk

	for bamfile in *sorted.bam
        	do
        	gatk --java-options "-Xmx12G" MarkDuplicates \
        	--INPUT $bamfile \
        	--OUTPUT $path/dedup_aligned_reads/$run/$(echo $bamfile | sed 's/sorted/duplicate_marked/') \
        	--METRICS_FILE $path/dedup_aligned_reads/$run/$(echo $bamfile | sed 's/bam/metrics/') \
        	--CREATE_INDEX true --ASSUME_SORTED true --TMP_DIR tmp --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
	done

module unload gatk

echo DEDUPLEXING FINISHED ON $(date) >> $path/log_files/$run_log.txt


## BASESCORE RECALIBRATION IN BAM FILES


reference='/disk/work/shared/data/genomes/ucsc.hg19_bwa_index/hg19.fa'
known1='/disk/work/shared/data/gatk/hg19_bgz/1000G_phase1.indels.b37_corrected.vcf'
known2='/disk/work/shared/data/gatk/hg19_bgz/Mills_and_1000G_gold_standard.indels.b37_corrected.vcf'
known3='/disk/work/shared/data/gatk/hg19_bgz/dbsnp_138.b37_corrected.vcf'

mkdir -p $path/recal_aligned_reads/$run
cd $path/dedup_aligned_reads/$run

module load gatk

	for bamfile in *.bam
        	do
		gatk BaseRecalibrator -I $bamfile -R $reference \
        	-O $path/recal_aligned_reads/$run/$(echo $bamfile | sed 's/duplicate_marked.bam/recal_data.table/') \
        	--known-sites $known1 \
        	--known-sites $known2 \
        	--known-sites $known3
		echo $(echo $bamfile | cut -f 1 -d "_") recalibration table finished >> $path/log_files/$run.log.txt
        done


echo RECALIBRATION TABLES FINISHED ON $(date) >> $path/log_files/$run.log.txt

	for bamfile in *.bam
        	do
		gatk ApplyBQSR -R $reference -I $bamfile \
        	--bqsr-recal-file $path/recal_aligned_reads/$run/$(echo $bamfile | sed 's/duplicate_marked.bam/recal_data.table/') \
        	-O $path/recal_aligned_reads/$run/$(echo $bamfile | sed 's/duplicate_marked.bam/recalibrated.bam/')
		echo $(echo $bamfile | cut -f 1 -d "_") recalibrated >> $path/log_files/$run.log.txt
        done

module unload gatk

echo FASTQ TO RECALIBRATED BAM FOR $run COMPLETED ON $(date) >> $path/log_files/$run.log.txt


reference='/disk/work/shared/data/genomes/ucsc.hg19_bwa_index/hg19.fa'
germline='/disk/work/shared/data/gnomad/af-only-gnomad.raw.sites_hg19.vcf'
somatic='/disk/work/shared/data/gatk/small_exac_common/small_exac_common_3.hg19.vcf'
pon='/disk/work/shared/data/gatk/pon/Mutect2-exome-panel.hg19.vcf'


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
reference='/disk/work/shared/data/genomes/ucsc.hg19_bwa_index/hg19.fa'

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
module load tabix

vcf2maf='/disk/work/shared/tools/vcf2maf/vcf2maf.pl'
reference='/disk/work/shared/data/genomes/ucsc.hg19_bwa_index/hg19.fa'
vep_exec='/disk/work/shared/tools/ensembl-vep'
vep_data='/disk/work/shared/data/vep'

cd $path/annotation/$run
mkdir -p $path/maf_files/$run

	for vcf in *vcf
        	do
		perl $vcf2maf --input-vcf $vcf \
        	--output-maf $path/maf_files/$run/$(echo $vcf | sed 's/annotated.vcf/maf/') \
        	--retain-info --retain-fmt \
        	--tumor-id $(echo $vcf | cut -f 1 -d "_") --inhibit-vep \
        	--ref-fasta $reference --vcf-tumor-id $(echo $vcf | cut -f 1 -d "_")

        	gzip $path/maf_files/$run/$(echo $vcf | sed 's/annotated.vcf/maf/')
		echo $(echo $bamfile | cut -f 1 -d "_") maf finished >> $path/log_files/$run.log.txt
        done

module unload samtools
module unload vcf2maf
module unload tabix

mkdir -p $path/tables/$run
cd $path/maf_files/$run

## CREATING VARIANT TABLES

Rscript='/disk/work/shared/tools/R-4.0.2/bin/Rscript'
xlsx='/disk/work/users/lv1/twist_myeloid_lymphoid/code/R/table_to_xlsx.R'
combined='/disk/work/users/lv1/twist_myeloid_lymphoid/code/R/combined_table.R'

for maffile in *.maf.gz
        do
        name=$(echo $maffile | cut -f 1 -d ".")

#        zcat $maffile > $path/tables/$run/$(echo $maffile | cut -f 1 -d ".").table
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

