#!/bin/bash

if [ ! $# == 1 ]; then
  echo "Usage: $0 <sample id>"
  echo "Or: parallel $0 ::: <sampleIDdir1> <sampleIDdir2> ... <sampleIDdirN>"
  exit
fi

TYPE=snps

# input <dir> parameter
PRJ=$1

# set working dir
WORK=~/workspace/vladimir/${PRJ}
if [ -d $WORK ]; then
  echo "${WORK}: Directory found."
else
  echo "${WORK}: Directory does not exist !"
  exit
fi

#set common prefix for processed files, e.g. /home/bla/bla .bam .vcf .log
PREFIX=${WORK}/${PRJ}

#set trio data file
PED=${PREFIX}.ped`
if [[ ! -f ${PED} ]]; then
	echo "${PED}: File not found."
	exit
fi

# get trio names from ped file 
F1=`awk 'NR == 1 { print $2 }' < ${PED}`
F2=`awk 'NR == 1 { print $3 }' < ${PED}`
F3=`awk 'NR == 1 { print $4 }' < ${PED}`
echo ${F1}
echo ${F2}
echo ${F3}


#
REF=ucsc.hg19
# or
#REF=hs37d5

RefPrefix=~/bin/${REF}

RefGENOME=${RefPrefix}/${REF}.fasta
# or
#RefGENOME=${RefPrefix}/${REF}.fa

#RefINTERVAL=${RefPrefix}/trusight_cancer_manifest.bed
# or
#RefINTERVAL=${RefPrefix}/target.intervals.bed

RefSNP=${RefPrefix}/dbsnp_137.${REF}.vcf
Ref1000G=${RefPrefix}/1000G_omni2.5.${REF}.vcf
RefHAPMAP=${RefPrefix}/hapmap_3.3.${REF}.vcf

#sw parameters:
#memory alloc.
MEM=-Xmx16g
NCORES=4
#bin
GATK=~/bin/gatk/GenomeAnalysisTK.jar
PSEQ=~/bin/pseq
#PICARD=~/bin/picard
#BWA=~/bin/bwa/bwa
#SAMTOOLS=~/bin/samtools/samtools

echo "$# parameters...";


#steps:
#
COMBINEALL=false
# or
PEDSNPCall=false

VQSR=false
# or
HARDFILTER=false

PHASETRANSM=false
RISKASSESS=false
RISKMODEL=false


if $COMBINEALL; then
	# merge all vcf
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T CombineVariants \
	--variant ${F1}.vcf \
	--variant ${F2}.vcf \
	--variant ${F3}.vcf \
	-o ${PREFIX}.all.vcf \
	-genotypeMergeOptions UNIQUIFY
fi


if $VQSR; then
	#8.1 quality score recalibration:
	# see also http://gatkforums.broadinstitute.org/discussion/2757/vqsr-for-small-exome-data-sets :
	#1)The 1000G samples should be included in the entire calling process, not just VQSR. 
	#This adds quite a bit of overhead of course, but it's necessary for the variant annotations to "fit".
	#2) The idea is to minimize differences by choosing samples from a similar cohort, e.g. match samples by geographical/ethnic origin and any other property possible.
	#It is also better if they were sequenced using a similar platform (so the error modes will be similar).

	# prepare snp model
	java ${MEM} -jar ${GATK} -T VariantRecalibrator \
		-nt $NCORES \
		-mode SNP \
		-R ${RefGENOME} \
		-input ${PREFIX}.all.vcf \
		--maxGaussians 4 \
		--resource:hapmap,known=false,	training=true,	truth=true,		prior=15.0	${RefHAPMAP} \
		--resource:omni,  known=false,	training=true,	truth=true,		prior=12.0	${Ref1000G} \
		--resource:dbsnp, known=true,		training=false,	truth=false,	prior=2.0	${RefSNP} \
		-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
		-recalFile ${PREFIX}.snps.recal \
		-tranchesFile ${PREFIX}.snps.tranches

	# apply SNP model
	java ${MEM} -jar ${GATK} -T ApplyRecalibration \
		-nt $NCORES \
		-mode SNP \
		-R ${RefGENOME} \
		-input ${PREFIX}.all.vcf \
		--ts_filter_level 99.0 \
		-tranchesFile ${PREFIX}.snps.tranches \
		-recalFile ${PREFIX}snps.recal \
		-o ${PREFIX}.snps.filt.vcf


	# prepare indel model
	java ${MEM} -jar ${GATK} -T VariantRecalibrator \
		-nt $NCORES \
		-mode INDEL \
		-R ${RefGENOME}
		-input ${PREFIX}.all.vcf \
		--maxGaussians 4 -std 10.0 -percentBad 0.12 \
		--resource:mills, known=true, training=true, truth=true, prior=12.0 ~/bin/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf \
		-an QD -an FS -an HaplotypeScore -an ReadPosRankSum \
		-recalFile ${PREFIX}.indels.recal \
		-tranchesFile ${PREFIX}.indels.tranches
	
	# apply indel model
	java ${MEM} -jar ${GATK} -T ApplyRecalibration \
		-nt $NCORES \
		-mode INDEL \
		-R ${RefGENOME} \
		-input ${PREFIX}.all.vcf \
		-ts_filter_level 99.9 \
		-recalFile ${PREFIX}.indels.recal \
		-tranchesFile ${PREFIX}.indels.tranches \
		-o ${PREFIX}.indels.filt.vcf
fi


if $HARDFILTER; then
	# Extract the SNPs from the call set
	java ${MEM} -jar ${GATK} -T SelectVariants \
		-nt $NCORES \
		-R ${RefGENOME} \
		-V ${PREFIX}.all.vcf \
		-selectType SNP \
		-o ${PREFIX}.snps.vcf
	#
	#
	# Apply the filter to the SNP call set
	java ${MEM} -jar ${GATK} -T VariantFiltration \
		-R ${RefGENOME} \
		-V ${PREFIX}.snps.vcf
		--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
		--filterName "my_snp_filter" \
		-o ${PREFIX}.snps.filt.vcf


	# Extract the Indels from the call set
	java ${MEM} -jar ${GATK} -T SelectVariants \
		-R ${RefGENOME} 
		-nt $NCORES \
		-V ${PREFIX}.all.vcf \
		-selectType INDEL \
		-o ${PREFIX}.indels.vcf
	#
	#
	# Apply the filter to the Indel call set
	java ${MEM} -jar ${GATK} -T VariantFiltration \
		-R ${RefGENOME} \
		-V ${PREFIX}.indels.vcf \
		--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
		--filterName "my_indel_filter" \
		-o ${PREFIX}.indels.filt.vcf
fi

if $PHASETRANSM; then
	#We recommend doing phasing after VQSR. (Geraldine: https://gatkforums.broadinstitute.org/discussion/2727/phasebytransmission-step-before-or-after-variantrecalibration)
	
	#Geraldine - http://gatkforums.broadinstitute.org/discussion/4192/about-trio-analysis:
	#That's right -- you run the preprocessing and calling step per sample, producing a gVCF for each member of the family, then you run joint genotyping on them together. No need to provide the ped file for calling or genotyping; at this time that information is not taken into account by the tools at those steps. You will use PhaseByTransmission (with pedigree file) after filtering.
	# Laurent - http://gatkforums.broadinstitute.org/discussion/47/phasebytransmission-with-more-than-just-trio:
	#Currently PhaseByTransmission only supports trios so you won't be able to use the information about all 4 siblings jointly.
	#Including just the trio of interest is the correct way to go for the moment, however if you leave the other siblings in the VCF file, you should either:
	# Add the children in the PED file but code them as unrelated individuals (they will simply be ignored)
	# Specify the flag --pedigreeValidationType SILENT. This flag lets the GATK run even if not all individuals are found in both the PED and VCF file.

	#snps trio/family
	java ${MEM} -jar ${GATK} -T PhaseByTransmission \
	-ped ${PED} \
	-V ${PREFIX}.snps.filt.vcf \
	-o ${PREFIX}.snps.phase.vcf
	#
	#snps quality metrics
	cd ~/bin/snpEff/
	cat ${PREFIX}.snps.phase.vcf | java -jar SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS' | FILTER =~ 'VQSRT')" > ${PREFIX}.snps.pass.vcf
	cd ${WORK}	
  ${PSEQ}/pseq ${PREFIX}.snps.phase.vcf v-stats > ${PREFIX}.snps.phase.stats
  ${PSEQ}/pseq ${PREFIX}.snps.pass.vcf  v-stats > ${PREFIX}.snps.pass.stats
  
  
  
	#indels trio
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T PhaseByTransmission 
	-ped ${PED} \
	-V ${PREFIX}.indels.filt.vcf \
	-o ${PREFIX}.indels.phase.vcf
	#
  #indels quality metrics
	cd ~/bin/snpEff/
	cat ${PREFIX}.indels.phase.vcf | java -jar SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS' | FILTER =~ 'VQSRT')" > ${PREFIX}.indels.pass.vcf
	cd ${WORK}	
	${PSEQ}/pseq ${PREFIX}.indels.phase.vcf v-stats > ${PREFIX}.indels.phase.stats 
	${PSEQ}/pseq ${PREFIX}.indels.pass.vcf  v-stats > ${PREFIX}.indels.pass.stats
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

	# snps section
	java ${MEM} -jar snpEff.jar \
	  -v \
	  -lof \
	  -motif \
	  -hgvs \
	  -nextProt \
	  GRCh37.74 \
	  ${PREFIX}.snps.pass.vcf > ${PREFIX}.snps.snpeff.vcf
	
	java ${MEM} -jar SnpSift.jar \
	  caseControl \
	  -v \
	  -tfam ${PED} \
	  ${PREFIX}.snps.snpeff.vcf > ${PREFIX}.snps.snpeff.cc.vcf

	# indels section
	java ${MEM} -jar snpEff.jar \
	  -v \
	  -lof \
	  -motif \
	  -hgvs \
	  -nextProt \
	  GRCh37.74 \
	  ${PREFIX}.indels.pass.vcf > ${PREFIX}.indels.snpeff.vcf
	
	java ${MEM} -jar SnpSift.jar \
	  caseControl \
	  -v \
	  -tfam ${PED} \
	  ${PREFIX}.indels.snpeff.vcf > ${PREFIX}.indels.snpeff.cc.vcf

fi

if $RISKMODEL; then
 	# filter hi moderate risk
 	#cat ${PREFIX}/${TYPE}/trio.snpeff_step2.vcf | java -jar SnpSift.jar filter "(EFF[*].IMPACT = 'HIGH')" 									> ${PREFIX}/${TYPE}/snpeff_risk_hi.vcf
	for TYPE in snps indels
	do
		cat ${PREFIX}.${TYPE}.snpeff.cc.vcf | java -jar SnpSift.jar filter "((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE'))" 	> ${PREFIX}.${TYPE}.snpeff.risk.vcf
		# rec. condition
		cat ${PREFIX}.${TYPE}.snpeff.risk.vcf | java -jar SnpSift.jar filter "(Cases[0] = 1) & (Controls[0] = 0)" 																								> ${PREFIX}.${TYPE}.snpeff.rec.vcf
		# dom. condition -> tutti de novo, se non ho un gen. malato
		cat ${PREFIX}.${TYPE}.snpeff.risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) & (Controls[1] = 0))" > ${PREFIX}.${TYPE}.snpeff.dom.vcf
		# mixed/polygenic model
		cat ${PREFIX}.${TYPE}.snpeff.risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) & ((Controls[1] = 2) | (Controls[1] = 1)))" > ${PREFIX}.${TYPE}.snpeff.hetmod.vcf
		#
		#split trio
		#vcf-subset -c ${F1}.variant1 ${WORK}.${TYPE}.pass.vcf > ${WORK}/${F1}.${TYPE}.pass.vcf
		#vcf-subset -c ${F2}.variant2 ${WORK}.${TYPE}.pass.vcf > ${WORK}/${F2}.${TYPE}.pass.vcf
		#vcf-subset -c ${F3}.variant3 ${WORK}.${TYPE}.pass.vcf > ${WORK}/${F3}.${TYPE}.pass.vcf
	done
fi


if false; then
	java -jar SnpSift.jar pedShow \
		-tfam ${PED} \
		${PREFIX}.${TYPE}.snpeff.risk.vcf \
		${WORK}/chart
	
	cat ${PREFIX}.${TYPE}.snpeff.risk.vcf | ./scripts/vcfInfoOnePerLine.pl
	
	java -Xmx1g -jar SnpSift.jar \
	    annotate \
	    -v \
	    protocols/db/clinvar_00-latest.vcf \
	    ${PREFIX}.${TYPE}.snpeff.risk.vcf > ${PREFIX}.${TYPE}.snpeff.risk.clinvar.vcf
	
	cat ${PREFIX}.${TYPE}.snpeff.risk.clinvar.vcf | ./scripts/vcfInfoOnePerLine.pl
fi