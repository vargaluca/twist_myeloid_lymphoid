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
		 /disk/work/shared/tools/fastqc/fastqc -o $path/raw_fastqc/$run $filename
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
		/disk/work/shared/tools/fastqc/fastqc -o $path/trimmed_fastqc/$run $filename
		echo $(echo $filename | cut -f 1 -d "_") trimmed fastqc finished >> $path/log_files/$run.log.txt
        done

cd $path/trimmed_fastqc/$run
multiqc . --filename $run-trimmed_multiqc
mv $run-trimmed_multiqc.html $path/tables/$run/

echo TRIMMED FASTQC FINISHED ON $(date) >> $path/log_files/$run.log.txt


## ALIGNING TRIMMED FASTQ TO REFERENCE GENOME WITH BWA-MEM


reference='/disk/work/shared/data/genomes/ucsc.hg19_bwa_index/ucsc.hg19.fasta.gz'
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

covered='/disk/work/diagnostics/TWLM/lv1/twist_all_aml_diagnostics/bed_files/Target_bases_covered_by_probes_Amplikon_LymphoidMyeloid_1X_TE-95491860_hg19_corrected.bed'

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


reference='/disk/work/shared/data/genomes/ucsc.hg19_bwa_index/ucsc.hg19.fasta.gz'
known1='/disk/work/shared/data/gatk/hg19_bgz/1000G_phase1.indels.hg19.sites.vcf.gz'
known2='/disk/work/shared/data/gatk/hg19_bgz/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz'
known3='/disk/work/shared/data/gatk/hg19_bgz/dbsnp_138.hg19.vcf.gz'

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

