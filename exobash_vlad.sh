#!/bin/bash

# align parameters:
# first input reads file
MAP1=$1
# second input reads file
MAP2=$2
# mapped output filename
MAPPED=$3
#
#
PRJ=cardio
TYPE=indels
WORK=~/workspace/vladimir
#RefGENOME=~/bin/hg19/ucsc.hg19.fasta
#RefPrefix=ucsc.hg19
RefGENOME=~/bin/hs37d5/hs37d5.fa
RefPrefix=hs37d5

#memory alloc.
MEM=-Xmx16g
NCORES=4
GATK=~/bin/gatk/GenomeAnalysisTK.jar
PSEQ=~/bin/pseq
PICARD=~/bin/picard

echo "$# parameters...";


#steps
QC=true
HGindex=false
MAP2Bam=false

PICARDPreproc=false
INDELRealign=false
#
# fixmate done by indelrealigner
FIXMATE=false
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


if $HGindex; then
	cd ~/bin/$RefPrefix
	# build index for bwa
	~/bin/bwa/bwa index -a bwtsw -p $RefPrefix $RefGENOME
	#mv ucsc.hg19* $HG19PATH
fi

if $QC; then
	#
	~/bin/FastQC/fastqc $(MAP1)
	~/bin/FastQC/fastqc $(MAP2)
	#
	# trimmomatic
fi

if $MAP2Bam; then
	# mapping su bundled-hg19 con bwa o bowtie
	#1 Align samples to genome (BWA), generates SAI files: .fastq ==> .sai
	~/bin/bwa/bwa aln -t $(NCORES) $(RefGENOME) -b $(MAP1).bam > $(MAP1).sai
	~/bin/bwa/bwa aln -t $(NCORES) $(RefGENOME) -b $(MAP2).bam > $(MAP2).sai
	
	#2 Convert SAI to SAM (BWA): r1.sai, r2.sai ==> .sam 
	~/bin/bwa/bwa sampe -P $(RefGENOME) $(MAP1).sai $(MAP2).sai $(MAP1).bam $(MAP2).bam > $(MAPPED).sam
	#this works with specific versions of bwa: 0.6.2/0.5.10 otherwise ==> picard/AddOrReplaceReadGroups
	#~/bin/bwa/bwa sampe -P $(RefGENOME) -r '@RG\tID:foo\tSM:bar' $(MAP1).sai $(MAP2).sai $(MAP1).bam $(MAP2).bam > $(MAPPED).sam
	
	#3 Convert SAM to BAM binary format (SAM Tools): .sam ==> .bam
	~/bin/samtools/samtools import $(RefGENOME).fai $(MAPPED).sam $(MAPPED).bam 

	#4 Sort BAM (SAM Tools, second param is a prefix): .bam ==> .sort.bam 
	~/bin/samtools/samtools sort $(MAPPED).bam $(MAPPED).sort

	#5 Index BAM (SAM Tools): .bam ==> .bai
	~/bin/samtools/samtools index $(MAPPED).sort.bam
fi

#
#   lexy.bash
#

if $PICARDPreproc; then
	#5.1 add grouping info in header (picard): .sort.bam ==> .group.bam, (.group.bai)
	java $MEM -Djava.io.tmpdir=/tmp -jar ${PICARD}/AddOrReplaceReadGroups.jar I=$(MAPPED).lexy.bam O=$(MAPPED).group.bam LB=whatever PL=illumina PU=r12 SM=$FASTQ1NAME CREATE_INDEX=true SO=coordinate VALIDATION_STRINGENCY=LENIENT

	#5.2 su cui marco i PCR duplicati group.bam ==> .marked.bam, .marked.bai
	java $MEM -Djava.io.tmpdir=/tmp -jar ${PICARD}/MarkDuplicates.jar INPUT=$(MAPPED).group.bam OUTPUT=$(MAPPED).marked.bam METRICS_FILE=$(MAPPED).metrics CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
	#??? create index x gatk ????
	#
	#java $MEM -Djava.io.tmpdir=/tmp -jar ${PICARD}/ReorderSam.jar I=$(MAPPED).marked.bam O=$(MAPPED).lexy_marked.bam REFERENCE= ${RefGENOME} CREATE_INDEX=true
fi

# here GATK needs .bai index
if $INDELRealign; then
	#6 Identify target regions for realignment (Genome Analysis Toolkit): .marked.bam ==> .intervals
	#  performance 500MB bam/~7000000 reads: 1 core -> 3hh 
	java $MEM -jar ${GATK} -T RealignerTargetCreator -nt $NCORES -R ${RefGENOME} -I $(MAPPED).marked.bam -o $(MAPPED).intervals
	
	#7 Realign BAM to get better Indel calling (Genome Analysis Toolkit): .marked.bam, .intervals ==> .realn.bam
	java -jar ${GATK} -T IndelRealigner -R ${RefGENOME} -I $(MAPPED).marked.bam -targetIntervals $(MAPPED).intervals -o $(MAPPED).realn.bam
fi

if $FIXMATE; then
	#7.1 using paired end data, the mate information must be fixed: .realn.bam ==> .mated.bam
	java $MEM -Djava.io.tmpdir=/tmp -jar ${PICARD}/FixMateInformation.jar INPUT=$(MAPPED).realn.bam OUTPUT=$(MAPPED).mated.bam SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
fi

if $MDTag; then
	#7.5 substantially improve SNP specificity with MD/NM editing: .mated.bam ==> .baq.bam
	~/bin/samtools/samtools calmd -Abr $(MAPPED).mated.bam ${RefGENOME} > $(MAPPED).baq.bam
	#8 Reindex the realigned BAM (SAM Tools)
	~/bin/samtools/samtools index $(MAPPED).baq.bam
fi

if $BASERecal; then
	#8.1 quality score recalibration: count covariates 
	java $MEM -jar ${GATK} \
	-l INFO \
	-R ${RefGENOME} \
	--DBSNP dbsnp_137.hg19.vcf \
	-I $(MAPPED).baq.bam \
	-T CountCovariates \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov DinucCovariate \
	-recalFile input.recal_data.csv

	#8.1 quality score recalibration: table recalibration
	java $MEM -jar ${GATK} \
	-l INFO \
	-R ${RefGENOME} \
	-I $(MAPPED).baq.bam \
	-T TableRecalibration \
	--out -I $(MAPPED).recal.bam \
	-recalFile input.recal_data.csv
fi

if $SNPCall; then
	#9 snp calls: .baq.bam ==> .vcf
	java $MEM -jar ${GATK} -glm BOTH -nt $NCORES -R ${RefGENOME} -T UnifiedGenotyper \
	-I $(MAPPED).baq.bam \
	-D ~/bin/hg19/dbsnp_137.hg19.vcf -metrics $(MAPPED).snps.metrics -A DepthOfCoverage -A AlleleBalance \
	-o $(MAPPED).snps.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 \
	-L ~/bin/hg19/target.intervals.bed
fi

if $PEDSNPCall; then
	#9 snp calls by family: .baq.bam ==> .vcf
	java $MEM -jar ${GATK} -glm BOTH -nt $NCORES $PED -R ${RefGENOME} -T UnifiedGenotyper \
	-I $(MAPPED).baq.bam \
	-D ~/bin/hg19/dbsnp_137.hg19.vcf -metrics $(MAPPED).snps.metrics -A DepthOfCoverage -A AlleleBalance \
	-o $(MAPPED).snps.vcf  -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 \
	-L ~/bin/hg19/target.intervals.bed
fi

if $COMBINEALL; then
	# merge all vcf
	java $MEM -jar ${GATK} -R ${RefGENOME} -T CombineVariants \
	--variant 5152F/C39G1ACXX_Pavlinic_lane5152F_sequence.snps.vcf \
	--variant 5388P/C39G1ACXX_Pavlinic_lane5388P_sequence.snps.vcf \
	--variant 5387M/C39G1ACXX_Pavlinic_lane5387M_sequence.snps.vcf \
	-o ${PRJ}.snps.vcf -genotypeMergeOptions UNIQUIFY
fi



if $VSQR; then
	#8.1 quality score recalibration: snp model
	java $MEM -jar ${GATK} -R ${RefGENOME} -T VariantRecalibrator -nt $NCORES \
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
	java $MEM -jar ${GATK} -R ${RefGENOME} -T VariantRecalibrator -nt $NCORES \
		-input $PRJ.snps.vcf \
		--maxGaussians 4 -std 10.0 -percentBad 0.12 \
		--resource:mills,known=true,training=true,truth=true,prior=12.0 ~/bin/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf \
		-an QD -an FS -an HaplotypeScore -an ReadPosRankSum \
		-mode INDEL \
		-recalFile $PRJ.indels.recal \
		-tranchesFile $PRJ.indels.tranches
	
	# apply SNP model
	java $MEM -jar ${GATK} -R ${RefGENOME} -T ApplyRecalibration -nt $NCORES \
		-input $PRJ.snps.vcf \
		--ts_filter_level 99.0 \
		-tranchesFile $PRJ.snps.tranches \
		-recalFile $PRJ.snps.recal \
		-mode SNP \
		-o $PRJ.snps.intermediate.vcf
	
	# apply indel model
	java $MEM -jar ${GATK} -R ${RefGENOME} -T ApplyRecalibration -nt $NCORES \
		-input $PRJ.snps.intermediate.vcf \
		-ts_filter_level 99.9 \
		-tranchesFile $PRJ.indels.tranches \
		-recalFile $PRJ.indels.recal \
		-mode INDEL \
		-o $PRJ.snps.ready.vcf
fi


if $HARDFILTER; then
	# Extract the SNPs from the call set
	java $MEM -jar ${GATK} -R ${RefGENOME} -T SelectVariants -nt $NCORES \
		-V $PRJ.snps.vcf \
	    -selectType SNP \
	    -o ${PRJ}/raw_snps.vcf
	#    
	#
	#
	# Apply the filter to the SNP call set
	java $MEM -jar ${GATK} -R ${RefGENOME} -T VariantFiltration \
	    -V ${PRJ}/raw_snps.vcf \
	    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
	    --filterName "my_snp_filter" \
	    -o ${PRJ}/filtered_snps.vcf

	# Extract the Indels from the call set
	java $MEM -jar ${GATK} -R ${RefGENOME} -T SelectVariants -nt $NCORES \
	-V $PRJ.snps.vcf \
	-selectType INDEL \
	-o ${PRJ}raw_indels.vcf
	#
	#
	# Apply the filter to the Indel call set
	java $MEM -jar ${GATK} -R ${RefGENOME} -T VariantFiltration \
	-V ${PRJ}/raw_indels.vcf \
	--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
	--filterName "my_indel_filter" \
	-o ${PRJ}/filtered_indels.vcf
fi

if $PHASETRANSM; then
	#snps trio
	java $MEM -jar ${GATK} -R ${RefGENOME} -T PhaseByTransmission -ped ${PRJ}/${PRJ}.ped \
	-V ${PRJ}/filtered_snps.vcf \
	-o ${PRJ}/trio.snps.phase.vcf

	# quality metrics
	cd ~/bin/snpEff/
	cat ${WORK}/${PRJ}/trio.snps.phase.vcf | java -jar SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS' | FILTER =~ 'VQSRT')" > ${WORK}/${PRJ}/trio.snps.pass.vcf
	cd ${WORK}	
  ${PSEQ}/pseq ${PRJ}//trio.snps.phase.vcf v-stats > ${PRJ}/snps.stats
  ${PSEQ}/pseq ${PRJ}//trio.snps.pass.vcf v-stats > ${PRJ}/snps.qc.stats
  
	#indels trio
	java $MEM -jar ${GATK} -R ${RefGENOME} -T PhaseByTransmission -ped ${PRJ}/${PRJ}.ped \
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
	#java $MEM -jar snpEff.jar eff -v -i vcf -o vcf hg19 ${WORK}/trio.snps.phase.vcf > ${WORK}/trio.snpeff.vcf
	# otherwise...
	# we annotate:
	# loss of function and nonsense mediated decay predictions (by adding the "-lof" command line option);
	# protein domain annotations from the curated NextProt database (option "-nextProt"); as well as
	# putative transcription factor binding sites from the ENSEMBL 'Regulatory Build' and Jaspar database (option "-motif"). 
	java $MEM -jar snpEff.jar \
	  -v \
	  -lof \
	  -motif \
	  -hgvs \
	  -nextProt \
	  GRCh37.74 \
	  ${WORK}/${PRJ}/trio.${TYPE}.pass.vcf > ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step1.vcf
	
	java $MEM -jar SnpSift.jar \
	  caseControl \
	  -v \
	  -tfam ${WORK}/${PRJ}/${PRJ}.ped \
	  ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step1.vcf > ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step2.vcf
fi

if $RISKMODEL; then
 	# filter hi moderate risk
	cat ${WORK}/${PRJ}/${TYPE}/trio.snpeff_step2.vcf | java -jar SnpSift.jar filter "((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE'))" 									> ${WORK}/${PRJ}/${TYPE}/snpeff_risk.vcf
  # rec. condition
	cat ${WORK}/${PRJ}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "(Cases[0] = 1) & (Controls[0] = 0)" 																								> ${WORK}/${PRJ}/${TYPE}/snpeff_rec.vcf
	# dom. condition
	#cat ${WORK}/${PRJ}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) | (Controls[1] = 0))" 			> ${WORK}/${PRJ}/${TYPE}/snpeff_dom.vcf
	cat ${WORK}/${PRJ}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) & (Controls[1] = 0))" 			> ${WORK}/${PRJ}/${TYPE}/snpeff_dom.vcf
fi

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
