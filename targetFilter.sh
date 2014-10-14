#!/bin/bash

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
RefGENOME=~/bin/hg19/ucsc.hg19.fasta
RefPrefix=ucsc.hg19
#RefGENOME=~/bin/hs37d5/hs37d5
#RefPrefix=hs37d5

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


#steps:
#
COMBINEALL=false
# or
PEDSNPCall=false

VSQR=false
# or
HARDFILTER=false

PHASETRANSM=false
RISKASSESS=false
RISKMODEL=false


if $COMBINEALL; then
	# merge all vcf
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T CombineVariants \
	--variant ${PREFIX}.vcf \
	--variant ${PREFIX}.vcf \
	--variant ${PREFIX}.vcf \
	-o ${PREFIX}.all.vcf \
	-genotypeMergeOptions UNIQUIFY
fi


if $VSQR; then
	#8.1 quality score recalibration:
	
	# prepare snp model
	java ${MEM} -jar ${GATK} -T VariantRecalibrator \
		-nt $NCORES \
		-mode SNP \
		-R ${RefGENOME} \
		-input ${PREFIX}.all.vcf \
		--maxGaussians 4 \
		--resource:hapmap,known=false,	training=true,	truth=true,		prior=15.0	~/bin/hg19/hapmap_3.3.hg19.vcf \
		--resource:omni,  known=false,	training=true,	truth=true,		prior=12.0	~/bin/hg19/1000G_omni2.5.hg19.vcf \
		--resource:dbsnp, known=true,		training=false,	truth=false,	prior=2.0		~/bin/hg19/dbsnp_137.hg19.vcf \
		-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
		-recalFile ${PREFIX}.snps.recal \
		-tranchesFile ${PREFIX}.snps.tranches

	# apply SNP model
	java ${MEM} -jar ${GATK} -T ApplyRecalibration \
		-nt $NCORES \
		-mode SNP \
		-R ${RefGENOME} \
		-input ${PREFIX}.vcf \
		--ts_filter_level 99.0 \
		-tranchesFile ${PREFIX}.snps.tranches \
		-recalFile ${PREFIX}snps.recal \
		-o ${PREFIX}.all.snps.vcf


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
		-input ${PREFIX}.snps.intermediate.vcf \
		-ts_filter_level 99.9 \
		-recalFile ${PREFIX}.indels.recal \
		-tranchesFile ${PREFIX}.indels.tranches \
		-o ${PREFIX}.all.indels.vcf
fi


if $HARDFILTER; then
	# Extract the SNPs from the call set
	java ${MEM} -jar ${GATK} -T SelectVariants \
		-nt $NCORES \
		-R ${RefGENOME} \
		-V ${PREFIX}.all.vcf \
		-selectType SNP \
		-o ${PREFIX}.raw.snps.vcf
	#
	#
	# Apply the filter to the SNP call set
	java ${MEM} -jar ${GATK} -T VariantFiltration \
		-R ${RefGENOME} \
		-V ${PREFIX}.raw.snps.vcf
		--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
		--filterName "my_snp_filter" \
		-o ${PREFIX}.all.snps.vcf


	# Extract the Indels from the call set
	java ${MEM} -jar ${GATK} -T SelectVariants \
		-R ${RefGENOME} 
		-nt $NCORES \
		-V ${PREFIX}.all.vcf \
		-selectType INDEL \
		-o ${PREFIX}.raw.indels.vcf
	#
	#
	# Apply the filter to the Indel call set
	java ${MEM} -jar ${GATK} -T VariantFiltration \
		-R ${RefGENOME} \
		-V ${PREFIX}.raw.indels.vcf \
		--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
		--filterName "my_indel_filter" \
		-o ${PREFIX}.all.indels.vcf
fi

if $PHASETRANSM; then
	#snps trio/family
	java ${MEM} -jar ${GATK} -T PhaseByTransmission \
		-ped ${PRJ}/${PRJ}.ped \
	-V ${PREFIX}.all.snps.vcf \
	-o ${PRJ}/trio.snps.phase.vcf

	# quality metrics
	cd ~/bin/snpEff/
	cat ${PREFIX}/trio.snps.phase.vcf | java -jar SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS' | FILTER =~ 'VQSRT')" > ${PREFIX}/trio.snps.pass.vcf
	cd ${WORK}	
  ${PSEQ}/pseq ${PRJ}//trio.snps.phase.vcf v-stats > ${PRJ}/snps.stats
  ${PSEQ}/pseq ${PRJ}//trio.snps.pass.vcf v-stats > ${PRJ}/snps.qc.stats
  
	#indels trio
	java ${MEM} -jar ${GATK} -R ${RefGENOME} -T PhaseByTransmission -ped ${PRJ}/${PRJ}.ped \
	-V ${PRJ}/filtered_indels.vcf \
	-o ${PRJ}/trio.indels.phase.vcf

  # quality metrics
	cd ~/bin/snpEff/
	cat ${PREFIX}/trio.indels.phase.vcf | java -jar SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS' | FILTER =~ 'VQSRT')" > ${PREFIX}/trio.indels.pass.vcf
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
	  ${PREFIX}/trio.${TYPE}.pass.vcf > ${PREFIX}/${TYPE}/trio.snpeff_step1.vcf
	
	java ${MEM} -jar SnpSift.jar \
	  caseControl \
	  -v \
	  -tfam ${PREFIX}/${PRJ}.ped \
	  ${PREFIX}/${TYPE}/trio.snpeff_step1.vcf > ${PREFIX}/${TYPE}/trio.snpeff_step2.vcf
fi

if $RISKMODEL; then
 	# filter hi moderate risk
 	#cat ${PREFIX}/${TYPE}/trio.snpeff_step2.vcf | java -jar SnpSift.jar filter "(EFF[*].IMPACT = 'HIGH')" 									> ${PREFIX}/${TYPE}/snpeff_risk_hi.vcf
	cat ${PREFIX}/${TYPE}/trio.snpeff_step2.vcf | java -jar SnpSift.jar filter "((EFF[*].IMPACT = 'HIGH') | (EFF[*].IMPACT = 'MODERATE'))" 									> ${PREFIX}/${TYPE}/snpeff_risk.vcf
  # rec. condition
	cat ${PREFIX}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "(Cases[0] = 1) & (Controls[0] = 0)" 																								> ${PREFIX}/${TYPE}/snpeff_rec.vcf
	# dom. condition -> tutti de novo, se non ho un gen. malato
	#cat ${PREFIX}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) | (Controls[1] = 0))" 			> ${PREFIX}/${TYPE}/snpeff_dom.vcf
	cat ${PREFIX}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) & (Controls[1] = 0))" 			> ${PREFIX}/${TYPE}/snpeff_dom.vcf
	# mixed/polygenic model
	cat ${PREFIX}/${TYPE}/snpeff_risk.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) | (Cases[1] = 1)) & ((Controls[0] = 0) & ((Controls[1] = 2) | (Controls[1] = 1)))" 			> ${PREFIX}/${TYPE}/snpeff_hetmod.vcf
fi
	#cat ${PREFIX}/${TYPE}/trio.snpeff_step2.vcf | java -jar SnpSift.jar filter "(EFF[*].IMPACT = 'HIGH')" 									> ${PREFIX}/${TYPE}/snpeff_risk_hi.vcf
	#cat ${PREFIX}/${TYPE}/snpeff_risk_hi.vcf | java -jar SnpSift.jar filter "((Cases[0] = 1) & (Controls[1] = 2))" 			> ${PREFIX}/${TYPE}/snpeff_hetmod.vcf
	#cat ${PREFIX}/${TYPE}/snpeff_hetmod.vcf | ./scripts/vcfInfoOnePerLine.pl
	
# split trio
#vcf-subset -c C39G1ACXX_Pavlinic_lane5152F_sequence.variant ${PREFIX}/trio.snps.pass.vcf > ${PREFIX}/5152F.snps.pass.vcf
#vcf-subset -c C39G1ACXX_Pavlinic_lane5388P_sequence.variant2 ${PREFIX}/trio.snps.pass.vcf > ${PREFIX}/5388P.snps.pass.vcf
#vcf-subset -c C39G1ACXX_Pavlinic_lane5387M_sequence.variant3 ${PREFIX}/trio.snps.pass.vcf > ${PREFIX}/5387M.snps.pass.vcf

#vcf-subset -c C39G1ACXX_Pavlinic_lane5152F_sequence.variant ${PREFIX}/trio.indels.pass.vcf > ${PREFIX}/5152F.indels.pass.vcf
#vcf-subset -c C39G1ACXX_Pavlinic_lane5388P_sequence.variant2 ${PREFIX}/trio.indels.pass.vcf > ${PREFIX}/5388P.indels.pass.vcf
#vcf-subset -c C39G1ACXX_Pavlinic_lane5387M_sequence.variant3 ${PREFIX}/trio.indels.pass.vcf > ${PREFIX}/5387M.indels.pass.vcf
	
if false; then	
	java -jar SnpSift.jar pedShow \
	  -tfam ${PREFIX}/${PRJ}.ped \
	  ${PREFIX}/${TYPE}/trio.snpeff_step3_ccriskhi.vcf \
	  ${WORK}/chart
	
	cat ${PREFIX}/${TYPE}/trio.snpeff_step3_ccriskhi.vcf | ./scripts/vcfInfoOnePerLine.pl
	
	java -Xmx1g -jar SnpSift.jar \
	    annotate \
	    -v \
	    protocols/db/clinvar_00-latest.vcf \
	    ${PREFIX}/${TYPE}/trio.snpeff_step3_ccriskhi.vcf > ${PREFIX}/${TYPE}/trio.snpeff_step3_ccriskhi.clinvar.vcf
	
		cat ${PREFIX}/${TYPE}/trio.snpeff_step3_ccriskhi.clinvar.vcf | ./scripts/vcfInfoOnePerLine.pl
fi

# back home
cd ${WORK}
