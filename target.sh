#!/bin/bash

# align parameters:
#
PRJ=$1
TYPE=snps
WORK=~/workspace/target/${PRJ}
# mapped output filename
FQ1=${WORK}/${PRJ}_R1.fastq
FQ2=${WORK}/${PRJ}_R2.fastq
MAP1=${WORK}/${PRJ}_R1.sai
MAP2=${WORK}/${PRJ}_R2.sai
MAPPED=${WORK}/${PRJ}.sam
#
RefGENOME=~/bin/hg19/ucsc.hg19.fasta
RefPrefix=ucsc.hg19
#RefGENOME=~/bin/hs37d5/hs37d5
#RefPrefix=hs37d5

#memory alloc.
MEM=-Xmx16g
NCORES=4
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
QC=false
HGindex=false
MAP2Bam=true
MAPFq2Bam=false
MAPBam2Bam=false

PICARDPreproc=false


INDELRealign=false
#
# FIXMATE step
# http://gatkforums.broadinstitute.org/discussion/1562/need-to-run-a-step-with-fixmateinformation-after-realignment-step
# Mark_DePristo Posts: 153 admin
# September 2012 Answer ?
# The indel realigner fixes mates itself, so the file is valid after realignment.
# How this is accomplished is a bit of magic code from Eric Banks and myself from several years ago.
# A previous version of the realigner did require you to run two passes, but that was before the magic code was written.
FIXMATE=false

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
MDTag=false

BASERecal=false
#
SNPCall=false

COMBINEALL=false
# or
PEDSNPCall=false

VSQR=false
# or
HARDFILTER=false

PHASETRANSM=false
RISKASSESS=false
RISKMODEL=false

#----------------------------------------------------------------------------------------------------------------------------
if $HGindex; then
	cd ~/bin/$RefPrefix
	# build index for bwa
	~/bin/bwa/bwa index -a bwtsw -p $RefPrefix $RefGENOME
	#mv ucsc.hg19* $HG19PATH
fi

#---------------------------------------------------------------------------------------------------------------------------
if $QC; then
	#
	~/bin/FastQC/fastqc --noextract ${MAP1}.fq ${MAP2}.fq
	#
	# trimmomatic
	java -jar ~/bin/trimmomatic/trimmomatic.jar PE -phred33 ${MAP1}.fq ${MAP2}.fq ${MAP1}_paired.fq ${MAP1}_unpaired.fq ${MAP2}_paired.fq ${MAP2}_unpaired.fq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	#
	~/bin/FastQC/fastqc --noextract 51960F_1_paired.fq 51960F_2_paired.fq
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
#  or                                                  -b .bam ==> .sai
if $MAP2BAM; then

	if $MAPBam2Bam; then
		#1.1 if reads in bam format...
		BAMFLAG=-b
	fi
	if $MAPFq2Bam; then
		#1.2 if reads in fastq format...
		BAMFLAG=
	fi

	${BWA} aln -t ${NCORES} ${RefGENOME} $BAMFLAG ${FQ1} > ${MAP1}
	${BWA} aln -t ${NCORES} ${RefGENOME} $BAMFLAG ${FQ2} > ${MAP2}
	#2 Convert SAI to SAM (BWA): r1.sai, r2.sai ==> .sam 
	${BWA} sampe -P ${RefGENOME} ${MAP1} ${MAP2} ${FQ1} ${FQ2} > ${MAPPED}
 	fi
	
	#this works with specific versions of bwa: 0.6.2/0.5.10 otherwise ==> picard/AddOrReplaceReadGroups
	#~/bin/bwa/bwa sampe -P ${RefGENOME} -r '@RG\tID:foo\tSM:bar' ${MAP1}.sai ${MAP2}.sai ${MAP1}.bam ${MAP2}.bam > ${MAPPED}.sam
	
	#3 Convert SAM to BAM binary format (SAM Tools): .sam ==> .bam
	${SAMTOOLS} import ${RefGENOME}.fai ${MAPPED} ${WORK}/${PRJ}.bam 

	#4 Sort BAM (SAM Tools, second param is a prefix): .bam ==> .sort.bam 
	${SAMTOOLS} sort ${WORK}/${PRJ}.bam ${WORK}/${PRJ}.sort

	#5 Index BAM (SAM Tools): .bam ==> .bai
	${SAMTOOLS} index ${WORK}/${PRJ}.sort.bam
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
	#5.1 add grouping info in header (picard): .sort.bam ==> .group.bam, (.group.bai)
	java ${MEM} -Djava.io.tmpdir=/tmp -jar ${PICARD}/AddOrReplaceReadGroups.jar I=${MAPPED}.bam O=${MAPPED}.group.bam LB=whatever PL=illumina PU=r12 SM=$FASTQ1NAME CREATE_INDEX=true SO=coordinate VALIDATION_STRINGENCY=LENIENT

	#5.2 su cui marco i PCR duplicati group.bam ==> .marked.bam, .marked.bai
	java ${MEM} -Djava.io.tmpdir=/tmp -jar ${PICARD}/MarkDuplicates.jar INPUT=${MAPPED}.group.bam OUTPUT=${MAPPED}.marked.bam METRICS_FILE=${MAPPED}.metrics CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
	#??? create index x gatk ????
	#
	#java ${MEM} -Djava.io.tmpdir=/tmp -jar ${PICARD}/ReorderSam.jar I=${MAPPED}.marked.bam O=${MAPPED}.lexy_marked.bam REFERENCE= ${RefGENOME} CREATE_INDEX=true
fi

# here GATK needs .bai index
if $INDELRealign; then
	#6 Identify target regions for realignment (Genome Analysis Toolkit): .marked.bam ==> .intervals
	#  performance 500MB bam/~7000000 reads: 1 core -> 3hh 
	java ${MEM} -jar ${GATK} -T RealignerTargetCreator -nt $NCORES \
		-R ${RefGENOME}.fa \
		-I ${MAPPED}.marked.bam \
		-o ${MAPPED}.intervals
		#--knonwn
		#-L
	
	#7 Realign BAM to get better Indel calling (Genome Analysis Toolkit): .marked.bam, .intervals ==> .realn.bam
	java -jar ${GATK} -T IndelRealigner -R ${RefGENOME}.fa -I ${MAPPED}.marked.bam -targetIntervals ${MAPPED}.intervals -o ${MAPPED}.realn.bam
fi

if $FIXMATE; then
	#7.1 using paired end data, the mate information must be fixed: .realn.bam ==> .mated.bam
	java ${MEM} -Djava.io.tmpdir=/tmp -jar ${PICARD}/FixMateInformation.jar INPUT=${MAPPED}.realn.bam OUTPUT=${MAPPED}.mated.bam SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
fi


if $MDTag; then
	#7.5 substantially improve SNP specificity with MD/NM editing: .mated.bam ==> .baq.bam
	~/bin/samtools/samtools calmd -Abr ${MAPPED}.mated.bam ${RefGENOME} > ${MAPPED}.baq.bam
	#8 Reindex the realigned BAM (SAM Tools)
	~/bin/samtools/samtools index ${MAPPED}.baq.bam
fi

if $BASERecal; then
	#8.1 quality score recalibration: count covariates 
	java ${MEM} -jar ${GATK} \
	-l INFO \
	-R ${RefGENOME} \
	--DBSNP dbsnp_137.hg19.vcf \
	-I ${MAPPED}.baq.bam \
	-T CountCovariates \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov DinucCovariate \
	-recalFile input.recal_data.csv

	#8.1 quality score recalibration: table recalibration
	java ${MEM} -jar ${GATK} \
	-l INFO \
	-R ${RefGENOME} \
	-I ${MAPPED}.baq.bam \
	-T TableRecalibration \
	--out -I ${MAPPED}.recal.bam \
	-recalFile input.recal_data.csv
fi


# lh3
# 08-27-2010, 01:19 PM
# Remember to set a quality threshold (about 50) on indels. Also the best indel caller so far is believed to be Dindel.
if $SNPCall; then
	#9 snp calls: .baq.bam ==> .vcf
	java ${MEM} -jar ${GATK} -glm BOTH -nt $NCORES -R ${RefGENOME} -T UnifiedGenotyper \
	-I ${MAPPED}.baq.bam \
	-D ~/bin/hg19/dbsnp_137.hg19.vcf -metrics ${MAPPED}.snps.metrics -A DepthOfCoverage -A AlleleBalance \
	-o ${MAPPED}.snps.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 \
	-L ~/bin/hg19/target.intervals.bed
fi

if $PEDSNPCall; then
	#9 snp calls by family: .baq.bam ==> .vcf
	java ${MEM} -jar ${GATK} -glm BOTH -nt $NCORES $PED -R ${RefGENOME} -T UnifiedGenotyper \
	-I ${MAPPED}.baq.bam \
	-D ~/bin/hg19/dbsnp_137.hg19.vcf -metrics ${MAPPED}.snps.metrics -A DepthOfCoverage -A AlleleBalance \
	-o ${MAPPED}.snps.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 \
	-L ~/bin/hg19/target.intervals.bed
fi

if $COMBINEALL; then
	# merge all vcf
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T CombineVariants \
	--variant 5152F/C39G1ACXX_Pavlinic_lane5152F_sequence.snps.vcf \
	--variant 5388P/C39G1ACXX_Pavlinic_lane5388P_sequence.snps.vcf \
	--variant 5387M/C39G1ACXX_Pavlinic_lane5387M_sequence.snps.vcf \
	-o ${PRJ}.snps.vcf -genotypeMergeOptions UNIQUIFY
fi



if $VSQR; then
	#8.1 quality score recalibration: snp model
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T VariantRecalibrator -nt $NCORES \
		-input $PRJ.snps.vcf \
		--maxGaussians 4 \
		--resource:hapmap,known=false,	training=true,	truth=true,		prior=15.0	~/bin/hg19/hapmap_3.3.hg19.vcf \
		--resource:omni,  known=false,	training=true,	truth=true,		prior=12.0	~/bin/hg19/1000G_omni2.5.hg19.vcf \
		--resource:dbsnp, known=true,		training=false,	truth=false,	prior=2.0		~/bin/hg19/dbsnp_137.hg19.vcf \
		-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
		-mode SNP \
		-recalFile $PRJ.snps.recal \
		-tranchesFile $PRJ.snps.tranches
	
	# indel model
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T VariantRecalibrator -nt $NCORES \
		-input $PRJ.snps.vcf \
		--maxGaussians 4 -std 10.0 -percentBad 0.12 \
		--resource:mills,known=true,training=true,truth=true,prior=12.0 ~/bin/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf \
		-an QD -an FS -an HaplotypeScore -an ReadPosRankSum \
		-mode INDEL \
		-recalFile $PRJ.indels.recal \
		-tranchesFile $PRJ.indels.tranches
	
	# apply SNP model
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T ApplyRecalibration -nt $NCORES \
		-input $PRJ.snps.vcf \
		--ts_filter_level 99.0 \
		-tranchesFile $PRJ.snps.tranches \
		-recalFile $PRJ.snps.recal \
		-mode SNP \
		-o $PRJ.snps.intermediate.vcf
	
	# apply indel model
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T ApplyRecalibration -nt $NCORES \
		-input $PRJ.snps.intermediate.vcf \
		-ts_filter_level 99.9 \
		-tranchesFile $PRJ.indels.tranches \
		-recalFile $PRJ.indels.recal \
		-mode INDEL \
		-o $PRJ.snps.ready.vcf
fi


if $HARDFILTER; then
	# Extract the SNPs from the call set
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T SelectVariants -nt $NCORES \
		-V $PRJ.snps.vcf \
	    -selectType SNP \
	    -o ${PRJ}/raw_snps.vcf
	#    
	#
	#
	# Apply the filter to the SNP call set
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T VariantFiltration \
	    -V ${PRJ}/raw_snps.vcf \
	    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
	    --filterName "my_snp_filter" \
	    -o ${PRJ}/filtered_snps.vcf

	# Extract the Indels from the call set
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T SelectVariants -nt $NCORES \
	-V $PRJ.snps.vcf \
	-selectType INDEL \
	-o ${PRJ}raw_indels.vcf
	#
	#
	# Apply the filter to the Indel call set
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T VariantFiltration \
	-V ${PRJ}/raw_indels.vcf \
	--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
	--filterName "my_indel_filter" \
	-o ${PRJ}/filtered_indels.vcf
fi

if $PHASETRANSM; then
	#snps trio
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T PhaseByTransmission -ped ${PRJ}/${PRJ}.ped \
	-V ${PRJ}/filtered_snps.vcf \
	-o ${PRJ}/trio.snps.phase.vcf

	# quality metrics
	cd ~/bin/snpEff/
	cat ${WORK}/${PRJ}/trio.snps.phase.vcf | java -jar SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS' | FILTER =~ 'VQSRT')" > ${WORK}/${PRJ}/trio.snps.pass.vcf
	cd ${WORK}	
  ${PSEQ}/pseq ${PRJ}//trio.snps.phase.vcf v-stats > ${PRJ}/snps.stats
  ${PSEQ}/pseq ${PRJ}//trio.snps.pass.vcf v-stats > ${PRJ}/snps.qc.stats
  
	#indels trio
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T PhaseByTransmission -ped ${PRJ}/${PRJ}.ped \
	-V ${PRJ}/filtered_indels.vcf \
	-o ${PRJ}/trio.indels.phase.vcf

  # quality metrics
	cd ~/bin/snpEff/
	cat ${WORK}/${PRJ}/trio.indels.phase.vcf | java -jar SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS' | FILTER =~ 'VQSRT')" > ${WORK}/${PRJ}/trio.indels.pass.vcf
	cd ${WORK}	
	${PSEQ}/pseq ${PRJ}/trio.indels.phase.vcf v-stats > ${PRJ}/indels.stats 
	${PSEQ}/pseq ${PRJ}/trio.indels.pass.vcf v-stats > ${PRJ}/indels.qc.stats
fi

cd ~/bin/snpEff/
if $RISKASSESS; then
	# from cureFFI pipeline
	#java ${MEM} -jar snpEff.jar eff -v -i vcf -o vcf hg19 ${WORK}/trio.snps.phase.vcf > ${WORK}/trio.snpeff.vcf
	# otherwise...
	# we annotate:
	# loss of function and nonsense mediated decay predictions (by adding the "-lof" command line option);
	# protein domain annotations from the curated NextProt database (option "-nextProt"); as well as
	# putative transcription factor binding sites from the ENSEMBL 'Regulatory Build' and Jaspar database (option "-motif"). 
	java ${MEM} -jar snpEff.jar \
	  -v \
	  -lof \
	  -motif \
	  -hgvs \
	  -nextProt \
	  GRCh37.74 \
	  ${WORK}/${PRJ}/trio.${TYPE}.pass.vcf > ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step1.vcf
	
	java ${MEM} -jar SnpSift.jar \
	  caseControl \
	  -v \
	  -tfam ${WORK}/${PRJ}/${PRJ}.ped \
	  ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step1.vcf > ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step2.vcf
fi

if $RISKMODEL; then
 	# filter hi moderate risk
 	#cat ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step2.vcf | java -jar SnpSift.jar filter "(EFF[*].IMPACT = 'HIGH')" 									> ${WORK}/${PRJ}/${TYPE}/snpeff_risk_hi.vcf
	cat ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step2.vcf | java -jar SnpSift.jar filter "((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE'))" 									> ${WORK}/${PRJ}/${TYPE}/snpeff_risk.vcf
  # rec. condition
	cat ${WORK}/${PRJ}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "(Cases[0] = 1) & (Controls[0] = 0)" 																								> ${WORK}/${PRJ}/${TYPE}/snpeff_rec.vcf
	# dom. condition -> tutti de novo, se non ho un gen. malato
	#cat ${WORK}/${PRJ}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) | (Controls[1] = 0))" 			> ${WORK}/${PRJ}/${TYPE}/snpeff_dom.vcf
	cat ${WORK}/${PRJ}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) & (Controls[1] = 0))" 			> ${WORK}/${PRJ}/${TYPE}/snpeff_dom.vcf
	# mixed/polygenic model
	cat ${WORK}/${PRJ}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) & ((Controls[1] = 2) | (Controls[1] = 1)))" 			> ${WORK}/${PRJ}/${TYPE}/snpeff_hetmod.vcf
fi
	#cat ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step2.vcf | java -jar SnpSift.jar filter "(EFF[*].IMPACT = 'HIGH')" 									> ${WORK}/${PRJ}/${TYPE}/snpeff_risk_hi.vcf
	#cat ${WORK}/${PRJ}/${TYPE}/snpeff_risk_hi.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) & (Controls[1] = 2))" 			> ${WORK}/${PRJ}/${TYPE}/snpeff_hetmod.vcf
	#cat ${WORK}/${PRJ}/${TYPE}/snpeff_hetmod.vcf | ./scripts/vcfInfoOnePerLine.pl
	
# split trio
#vcf-subset -c C39G1ACXX_Pavlinic_lane5152F_sequence.variant ${WORK}/${PRJ}/trio.snps.pass.vcf > ${WORK}/${PRJ}/5152F.snps.pass.vcf
#vcf-subset -c C39G1ACXX_Pavlinic_lane5388P_sequence.variant2 ${WORK}/${PRJ}/trio.snps.pass.vcf > ${WORK}/${PRJ}/5388P.snps.pass.vcf
#vcf-subset -c C39G1ACXX_Pavlinic_lane5387M_sequence.variant3 ${WORK}/${PRJ}/trio.snps.pass.vcf > ${WORK}/${PRJ}/5387M.snps.pass.vcf

vcf-subset -c C39G1ACXX_Pavlinic_lane5152F_sequence.variant ${WORK}/${PRJ}/trio.indels.pass.vcf > ${WORK}/${PRJ}/5152F.indels.pass.vcf
vcf-subset -c C39G1ACXX_Pavlinic_lane5388P_sequence.variant2 ${WORK}/${PRJ}/trio.indels.pass.vcf > ${WORK}/${PRJ}/5388P.indels.pass.vcf
vcf-subset -c C39G1ACXX_Pavlinic_lane5387M_sequence.variant3 ${WORK}/${PRJ}/trio.indels.pass.vcf > ${WORK}/${PRJ}/5387M.indels.pass.vcf
	
if false; then	
	java -jar SnpSift.jar pedShow \
	  -tfam ${WORK}/${PRJ}/${PRJ}.ped \
	  ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step3_ccriskhi.vcf \
	  ${WORK}/chart
	
	cat ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step3_ccriskhi.vcf | ./scripts/vcfInfoOnePerLine.pl
	
	java -Xmx1g -jar SnpSift.jar \
	    annotate \
	    -v \
	    protocols/db/clinvar_00-latest.vcf \
	    ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step3_ccriskhi.vcf > ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step3_ccriskhi.clinvar.vcf
	
		cat ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step3_ccriskhi.clinvar.vcf | ./scripts/vcfInfoOnePerLine.pl
fi

# back home
cd ${WORK}
