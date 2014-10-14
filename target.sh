#!/bin/bash
# e.g.
# ./target.sh ::: <sampleIDdir>
# or
# parallel ./target.sh ::: <sampleIDdir1> <sampleIDdir2> ... <sampleIDdirN>

if [ ! $# == 1 ]; then
  echo "Usage: $0 <sample id>"
  echo "Or: parallel $0 ::: <sampleIDdir1> <sampleIDdir2> ... <sampleIDdirN>"
  exit
fi

weight="$1"

TYPE=snps

PRJ=$1

WORK=~/workspace/target/${PRJ}
if [ -d $WORK ]; then
  echo "${WORK}: Directory found."
else
  echo "${WORK}: Directory does not exist !"
  exit
fi


PREFIX=${WORK}/${PRJ}
# mapped output filename
FQ1=${PREFIX}_R1.fastq
FQ2=${PREFIX}_R2.fastq
#
#
REF=ucsc.hg19
# or
#REF=hs37d5

RefPrefix=~/bin/${REF}

RefGENOME=${RefPrefix}/${REF}.fasta
# or
#RefGENOME=${RefPrefix}/${REF}.fa

RefINTERVAL=${RefPrefix}/trusight_cancer_manifest.bed
# or
RefINTERVAL=${RefPrefix}/target.intervals.bed

RefSNP=${RefPrefix}/dbsnp_137.${REF}.vcf

#align parameters:
#memory alloc.
MEM=-Xmx16g
NCORES=4
#bin
GATK=~/bin/gatk/GenomeAnalysisTK.jar
PSEQ=~/bin/pseq
PICARD=~/bin/picard
BWA=~/bin/bwa/bwa
SAMTOOLS=~/bin/samtools/samtools

echo "$# parameters...";


# general comments:
# lh3
# 10-17-2010, 08:31 PM
# The following is what I would recommend for Illumina:
# 1. Do alignment with novoalign or bwa. Mosaik is also great, but unfortunately it does not write soft clippings, which will affect programs at a later step. Bowtie is not recommended because it does not do gapped alignment.
# 2. If you have bandwidth, do realignment with GATK. If you do not, it actually does not matter too much. The major downside of not doing realignment is you may get confusing alignment in an alignment viewer (the most frequent question is "why the indel caller is calling an indel when there is only 1 read supporting that?").
# 3. Cap base quality BAQ (with samtools).
# 4. Call SNPs with whatever SNP caller. It does not matter too much when indels are cleaned.
# 5. Call indels with dindel. The Dindel group has shown convincing evidence that it is clearly better. When I evaluate it by myself, I am convinced again. It is much more sensitive than gatk realigner+IndelGenotyperV2; its specificity is also better. For exonome sequencing the difference is probably smaller because the hard regions are mostly related to repeats.


# http://www.biostars.org/p/1268/
# Question: What is the best pipeline for human whole exome sequencing?


#steps

HGindex=false
#
QC=false

# if input reads in bam format:
# BAMFLAG=-b

# or
MAP2Bam=true
#
PICARDPreproc=true
#
INDELRealign=true
#
# FIXMATE step
# http://gatkforums.broadinstitute.org/discussion/1562/need-to-run-a-step-with-fixmateinformation-after-realignment-step
# Mark_DePristo Posts: 153 admin
# September 2012 Answer ?
# The indel realigner fixes mates itself, so the file is valid after realignment.
# How this is accomplished is a bit of magic code from Eric Banks and myself from several years ago.
# A previous version of the realigner did require you to run two passes, but that was before the magic code was written.
# => FIXMATE=false

# Michael.James.Clark
# 10-26-2010, 11:41 AM
# Hi Heng,
# 3. Cap base quality BAQ (with samtools).
# What do you mean by "cap base quality"?
# Can you give more detailed suggestions about how to effectively use BAQ?
# We would like to include it in our pipeline but are unsure about how best to utilize it.
#
# drio
# 10-27-2010, 04:53 AM
# Here. http://seqanswers.com/forums/showpost.php?p=27334&postcount=26
# and here http://sourceforge.net/p/samtools/mailman/samtools-help/thread/EEEA3872-183F-4C2B-8FFB-7CE3EE877303@sanger.ac.uk/
# MDTag=false
FIXMATE_AND_BAQ=false
# or
BASERecal=true
#
SNPCall=true


#----------------------------------------------------------------------------------------------------------------------------
if $HGindex; then
	#cd ${RefPrefix}
	# build index for bwa
	${BWA} index -a bwtsw -p $RefPrefix $RefGENOME
	#mv ucsc.hg19* $HG19PATH
fi

#---------------------------------------------------------------------------------------------------------------------------
if $QC; then
	#
	if [[ ! -f ${FQ1} ]]; then
		echo "${FQ1} does not exist !"
	  exit
	fi
	if [[ ! -f ${FQ2} ]]; then
		echo "${FQ2} does not exist !"
	  exit
	fi

	~/bin/FastQC/fastqc --noextract ${FQ1} ${FQ2}
	#
	# trimmomatic
	java -jar ~/bin/trimmomatic/trimmomatic.jar PE -phred33 ${FQ1} ${FQ2} ${PREFIX}_R1_paired.fq ${PREFIX}_R1_unpaired.fq ${PREFIX}_R2_paired.fq ${PREFIX}_R2_unpaired.fq ILLUMINACLIP:~/bin/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:150
	#
	~/bin/FastQC/fastqc --noextract ${PREFIX}_R1_paired.fq ${PREFIX}_R2_paired.fq
	FQ1=${PREFIX}_R1_paired.fq
	FQ2=${PREFIX}_R2_paired.fq
	echo "Fastq reads updated to <${FQ1},${FQ2}>"
fi


# further justification for BWA mem
# http://www.bioplanet.com/gcat/reports/80/alignment/250bp-pe-small-indel/bwa/compare-68-72

# http://lists.idyll.org/pipermail/ngs-2013/2013-June/000099.html
#On Tue, Jun 25, 2013 at 2:04 PM, Rayan Chikhi <rchikhi at gmail.com> wrote:
#
# Yes, the aln/sampe commands invoke the original BWA algorithm
# (backtracking, from 2009), and the mem command invokes the latest BWA
# improvement (MEM algorithm, 2013).
# See http://bio-bwa.sourceforge.net/bwa.shtml#2
# I think it is fair to say that mem should replace aln/sampe for everything
# now. (except maybe on older, < 70 bp reads)
#
# Rayan
# mapping su bundled-hg19 con bwa o bowtie
#1 Align samples to genome (BWA), generates SAI files: .fastq ==> .sai
if $MAP2Bam; then
	if [[ ! -f ${FQ1} ]]; then
		echo "${FQ1} does not exist !"
	  exit
	fi
	if [[ ! -f ${FQ2} ]]; then
		echo "${FQ2} does not exist !"
	  exit
	fi

	${BWA} aln -t ${NCORES} ${RefGENOME} ${FQ1} > ${PREFIX}_R1.sai
	${BWA} aln -t ${NCORES} ${RefGENOME} ${FQ2} > ${PREFIX}_R2.sai
	
#2 Convert SAI to SAM (BWA): r1.sai, r2.sai ==> .sam 
	${BWA} sampe -P ${RefGENOME} ${PREFIX}_R1.sai ${PREFIX}_R2.sai ${FQ1} ${FQ2} > ${PREFIX}.sam
	
	# LH3:
	#this works with specific versions of bwa: 0.6.2/0.5.10 otherwise ==> picard/AddOrReplaceReadGroups
	#~/bin/bwa/bwa sampe -P ${RefGENOME} -r '@RG\tID:foo\tSM:bar' ${MAP1}.sai ${MAP2}.sai ${MAP1}.bam ${MAP2}.bam > ${PREFIX}.sam
	
	#3 Convert SAM to BAM binary format (SAM Tools): .sam ==> .bam
	${SAMTOOLS} import ${RefGENOME}.fai ${PREFIX}.sam ${PREFIX}.bam 

	#4 Sort BAM (SAM Tools, second param is a prefix): .bam ==> .sort.bam 
	${SAMTOOLS} sort ${PREFIX}.bam ${PREFIX}.sort

	#5 Index BAM (SAM Tools): .bam ==> .bai
	${SAMTOOLS} index ${PREFIX}.sort.bam
fi

#
#   lexy.bash
#

# http://www.broadinstitute.org/gatk/guide/best-practices#tutorials_mapdedup2909
#Fix a BAM that is not indexed or not sorted, has not had duplicates marked, or is lacking read group information.
#These steps can be performed independently of each other but this order is recommended.
#Prerequisites:  Installed Picard tools
#
#Steps:
#    Sort the aligned reads by coordinate order
#    Mark duplicates
#    Add read group information
#    Index the BAM file

if $PICARDPreproc; then
	#5.1 add grouping info in header (picard): .sort.bam ==> .group.bam, (.group.bai ?)
	java ${MEM} -Djava.io.tmpdir=/tmp -jar ${PICARD}/AddOrReplaceReadGroups.jar \
		I=${PREFIX}.sort.bam \
		O=${PREFIX}.group.bam \
		LB=whatever PL=illumina PU=r12 SM=$FASTQ1NAME SO=coordinate \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=LENIENT
# create index forse inutile qui ==> meglio dopo (5.2)

	#5.2 su cui marco i PCR duplicati group.bam ==> .marked.bam, .marked.bai
	java ${MEM} -Djava.io.tmpdir=/tmp -jar ${PICARD}/MarkDuplicates.jar \
		INPUT=${PREFIX}.group.bam \
		OUTPUT=${PREFIX}.marked.bam \
		METRICS_FILE=${PREFIX}.metrics \
		ASSUME_SORTED=true \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=LENIENT
	#??? create index x gatk ????
	#
	#java ${MEM} -Djava.io.tmpdir=/tmp -jar ${PICARD}/ReorderSam.jar I=${PREFIX}.marked.bam O=${PREFIX}.lexy_marked.bam REFERENCE= ${RefGENOME} CREATE_INDEX=true
fi

# here GATK needs .bai index
if $INDELRealign; then
	#6 Identify target regions for realignment (Genome Analysis Toolkit): .marked.bam ==> .intervals
	#  performance 500MB bam/~7000000 reads: 1 core -> 3hh 
	java ${MEM} -jar ${GATK} -T RealignerTargetCreator \
		-nt $NCORES \
		-R ${RefGENOME} \
		-I ${PREFIX}.marked.bam \
		-o ${PREFIX}.intervals
		#--knonwn
		#-L
	
	#7 Realign BAM to get better Indel calling (Genome Analysis Toolkit): .marked.bam, .intervals ==> .realn.bam
	java -jar ${GATK} -T IndelRealigner \
		-R ${RefGENOME} \
		-I ${PREFIX}.marked.bam \
		-targetIntervals ${PREFIX}.intervals \
		-o ${PREFIX}.realn.bam
		# -baq CALCULATE_AS_NECESSARY
fi

# Geraldine (here: http://gatkforums.broadinstitute.org/discussion/1645/is-picard-fixmate-necessary) said:
# G:"The realigner handles this automatically for you, no need to run fixmates."
# -"Hello, is then necessary fixmate necessary in any step of the standard gatk protocol?"
# G:"No, normally you shouldn't need to use fixmates."
if $FIXMATE_AND_BAQ; then
	#7.1 using paired end data, the mate information must be fixed: .realn.bam ==> .mated.bam
	java ${MEM} -Djava.io.tmpdir=/tmp -jar ${PICARD}/FixMateInformation.jar \
		INPUT=${PREFIX}.realn.bam \
		OUTPUT=${PREFIX}.mated.bam \
		SO=coordinate \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=LENIENT

	#7.2 if using samtools => disable GATK baq
	GATK_BAQ_POLICY=OFF
	
	#7.3 substantially improve SNP specificity with MD/NM editing: .mated.bam ==> .baq.bam
	${SAMTOOLS} calmd -Abr ${PREFIX}.mated.bam ${RefGENOME} > ${PREFIX}.baq.bam
	#7.4 Reindex the realigned BAM (SAM Tools)
	${SAMTOOLS} index ${PREFIX}.baq.bam
else
	# trust the GATK pipeline for mate and baq
	GATK_BAQ_POLICY=RECALCULATE
fi

if $BASERecal; then
	#8.1 base quality score recalibration: count covariates 
	java ${MEM} -jar ${GATK} -T BaseRecalibrator \
		-R ${RefGENOME} \
		--knownSites ${RefSNP} \
		-cov ReadGroupCovariate \
		-cov QualityScoreCovariate \
		-cov CycleCovariate \
		-I ${PREFIX}.realn.bam \
		-o ${PREFIX}.recal_data.table
#		retired option:
#		-cov DinucCovariate \

# 8.2 output recalibrated bam
# Prints the first 2000 reads in the BAM file
	java ${MEM} -jar ${GATK} -T PrintReads \
		-R ${RefGENOME} \
		--BQSR ${PREFIX}.recal_data.table\
		-I ${PREFIX}.realn.bam \
		-o ${PREFIX}.baq.bam
# riaggiungere con GATK > 2.4
# vedi http://gatkforums.broadinstitute.org/discussion/2267/baq-tag-error
#		-baq ${GATK_BAQ_POLICY} \ 
fi 


# lh3
# 08-27-2010, 01:19 PM
# Remember to set a quality threshold (about 50) on indels. Also the best indel caller so far is believed to be Dindel.
if $SNPCall; then
	#9 snp calls: .baq.bam ==> .vcf
	java ${MEM} -jar ${GATK} -T UnifiedGenotyper \
		-glm BOTH \
		-nt $NCORES \
		-R ${RefGENOME} \
		-D ${RefSNP} \
		-metrics ${PREFIX}.snps.metrics \
		-A DepthOfCoverage \
		-A AlleleBalance \
		-stand_call_conf 50.0 \
		-stand_emit_conf 10.0 \
		-dcov 1000 \
		-L ${RefINTERVAL} \
		-baq CALCULATE_AS_NECESSARY \
		-I ${PREFIX}.baq.bam \
		-o ${PREFIX}.vcf
fi
