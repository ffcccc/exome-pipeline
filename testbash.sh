#!/bin/bash
FASTQPATH=$1
FASTQ1NAME=$2
FASTQ2NAME=$3
echo "$# parameters"; 

BIN="~/bin/"
echo
echo "BIN path: "$BIN
HG19PATH=$BIN"hg19/"
echo
echo "HG19 path: "$HG19PATH
HG19fa=$HG19PATH"ucsc.hg19.fasta" 
#HG19fa=$HG19PATH"hg19.fa" 
echo "reference genome: "$HG19fa
echo
echo "working dir: "$FASTQPATH
TRIM1fq=$FASTQPATH$FASTQ1NAME
echo
echo "working file1: "$TRIM1fq"("$(wc -l $TRIM1fq)"reads)"
TRIM2fq=$FASTQPATH$FASTQ2NAME
echo
echo "working file2: "$TRIM2fq"("$(wc -l $TRIM2fq)"reads)"

#memory alloc.
MEM=-Xmx16g

#steps
QC=false
HGindex=false
MAP2Bam=true
PICARDPreproc=true
INDELRealign=true
MDTag=true
SCORERecal=false
SNPCall=true


if $QC; then
	#
	~/bin/FastQC/fastqc $FASTQ1NAME &
	~/bin/FastQC/fastqc $FASTQ2NAME &
	#
	~/bin/fastx/bin/fastx_trimmer -Q33 -f 1 -l 95 -i $FASTQ1NAME -o $TRIM1 &
	~/bin/fastx/bin/fastx_trimmer -Q33 -f 1 -l 95 -i $FASTQ2NAME -o $TRIM2 &

	#check for bad qc reads and exclude them...
	#Illumina reads are often pre-processed to remove reads based on some quality filtering stratergy.
	#For example using the fastq_quality_filter tool of the FASTX-Toolkit to remove reads that do not have a minimum quality value of 20 across the whole read:
	#e.g.
	#fastq_quality_filter  -q 20 –p 100 -i s_1_1_sequence.txt -o s_1_1_sequence.txt_filtered_q20_p100.fastq
	#fastq_quality_filter  -q 20 –p 100 -i s_1_2_sequence.txt -o s_1_2_sequence.txt_filtered_q20_p100.fastq
	#
	#This type of filtering can and will result in reads that are nolonger paired.

	~/bin/fastx/bin/fastq_quality_filter -Q33 -q 20 –p 100 -i $TRIM1.fastq -o 111337/111337_R1_trim_q20.fastq &
	~/bin/fastx/bin/fastq_quality_filter -Q33 -q 20 –p 100 -i $TRIM2.fastq -o 111337/111337_R2_trim_q20.fastq &

	~/bin/FastQC/fastqc 111337/111337_R1_trim_q20.fastq &
	~/bin/FastQC/fastqc 111337/111337_R2_trim_q20.fastq &
fi

if $HGindex; then
	# build index for bwa
	~/bin/bwa/bwa index -a bwtsw -p ucsc.hg19 $HG19fa
	#mv ucsc.hg19* $HG19PATH
fi

if $MAP2Bam; then
	# mapping su bundled-hg19 con bwa o bowtie
	#1 Align samples to genome (BWA), generates SAI files: .fastq ==> .sai
	~/bin/bwa/bwa aln -t 4 ~/bin/hg19/ucsc.hg19 $TRIM1fq > $FASTQPATH"MAP1.sai"
	~/bin/bwa/bwa aln -t 4 ~/bin/hg19/ucsc.hg19 $TRIM2fq > $FASTQPATH"MAP2.sai"
	
	#2 Convert SAI to SAM (BWA): r1.sai, r2.sai ==> .sam 
	~/bin/bwa/bwa sampe -P ~/bin/hg19/ucsc.hg19 $FASTQPATH"MAP1.sai" $FASTQPATH"MAP2.sai" $TRIM1fq $TRIM2fq > $FASTQPATH"MAP.sam"	
	#this works with specific versions of bwa: 0.6.2/0.5.10 otherwise ==> picard/AddOrReplaceReadGroups
	#~/bin/bwa/bwa sampe -P ~/bin/hg19/ucsc.hg19 -r '@RG\tID:foo\tSM:bar' $FASTQPATH"MAP1.sai" $FASTQPATH"MAP2.sai" $TRIM1fq $TRIM2fq > $FASTQPATH"MAP.sam"
	
	#3 Convert SAM to BAM binary format (SAM Tools): .sam ==> .bam
	~/bin/samtools/samtools import $HG19PATH"ucsc.hg19.fasta.fai" $FASTQPATH"MAP.sam" $FASTQPATH"MAP.bam" 

	#4 Sort BAM (SAM Tools, second param is a prefix): .bam ==> .sort.bam 
	~/bin/samtools/samtools sort $FASTQPATH"MAP.bam" $FASTQPATH"MAP.sort"

	#5 Index BAM (SAM Tools): .bam ==> .bai
	~/bin/samtools/samtools index $FASTQPATH"MAP.sort.bam"

if $PICARDPreproc; then
	#5.1 add grouping info in header (picard): .sort.bam ==> .group.bam, .group.bai
	java $MEM -Djava.io.tmpdir=/tmp -jar ~/bin/picard/AddOrReplaceReadGroups.jar I=$FASTQPATH"MAP.sort.bam" O=$FASTQPATH"MAP.group.bam" LB=whatever PL=illumina PU=r12 SM=111337 CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

	#5.2 su cui marco i PCR duplicati group.bam ==> .marked.bam, .marked.bai
	java $MEM -Djava.io.tmpdir=/tmp -jar ~/bin/picard/MarkDuplicates.jar INPUT=$FASTQPATH"MAP.group.bam" OUTPUT=$FASTQPATH"MAP.marked.bam" METRICS_FILE=metrics CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
fi

if $INDELRealign; then
	#6 Identify target regions for realignment (Genome Analysis Toolkit): .marked.bam ==> .intervals 
	java $MEM -jar ~/bin/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/bin/hg19/ucsc.hg19.fasta -I $FASTQPATH"MAP.marked.bam" -o $FASTQPATH"MAP.intervals"
	
	#7 Realign BAM to get better Indel calling (Genome Analysis Toolkit): .marked.bam, .intervals ==> .realn.bam
	java -jar ~/bin/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R ~/bin/hg19/ucsc.hg19.fasta -I $FASTQPATH"MAP.marked.bam" -targetIntervals $FASTQPATH"MAP.intervals" -o $FASTQPATH"MAP.realn.bam"
	
	#7.1 using paired end data, the mate information must be fixed: .realn.bam ==> .mated.bam
	java $MEM -Djava.io.tmpdir=/tmp -jar ~/bin/picard/FixMateInformation.jar INPUT=$FASTQPATH"MAP.realn.bam" OUTPUT=$FASTQPATH"MAP.mated.bam" SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
fi

if $MDTag; then
	#7.5 substantially improve SNP specificity with MD/NM editing: .mated.bam ==> .baq.bam"
	~/bin/samtools/samtools calmd -Abr $FASTQPATH"MAP.mated.bam" ~/bin/hg19/ucsc.hg19.fasta > $FASTQPATH"MAP.baq.bam"
	#8 Reindex the realigned BAM (SAM Tools)
	~/bin/samtools/samtools index $FASTQPATH"MAP.baq.bam"
fi

if $SCORERecal; then
	#8.1 quality score recalibration: count covariates 
	java $MEM -jar ~/bin/gatk/GenomeAnalysisTK.jar \
	-l INFO \
	-R ~/bin/hg19/ucsc.hg19.fasta \
	--DBSNP dbsnp_137.hg19.vcf \
	-I $FASTQPATH"MAP.baq.bam" \
	-T CountCovariates \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov DinucCovariate \
	-recalFile input.recal_data.csv

	#8.1 quality score recalibration: table recalibration
	java $MEM -jar ~/bin/gatk/GenomeAnalysisTK.jar \
	-l INFO \
	-R ~/bin/hg19/ucsc.hg19.fasta \
	-I $FASTQPATH"MAP.baq.bam" \
	-T TableRecalibration \
	--out -I $FASTQPATH"MAP.recal.bam" \
	-recalFile input.recal_data.csv
fi

if $SNPCall; then
	#9 snp calls: .baq.bam ==> .vcf
	java $MEM -jar ~/bin/gatk/GenomeAnalysisTK.jar -glm BOTH -R ~/bin/hg19/ucsc.hg19.fasta -T UnifiedGenotyper -I $FASTQPATH"MAP.baq.bam" -D ~/bin/hg19/dbsnp_137.hg19.vcf -metrics $FASTQPATH"MAP.snps.metrics" -A DepthOfCoverage -A AlleleBalance -o $FASTQPATH"MAP.snps.vcf"  -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 -L ~/bin/hg19/target.intervals.bed
fi

