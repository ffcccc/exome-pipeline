#!/bin/bash
BAMPATH=$1
BAMNAME=$2
MEM=-Xmx16g

~/bin/samtools/samtools view -H $BAMPATH"/"$BAMNAME".sort.bam" > $BAMPATH"/header.sam"
sed "s/.fa/ /" $BAMPATH"/header.sam" > $BAMPATH"/header_corrected.sam"
~/bin/samtools/samtools reheader $BAMPATH"/header_corrected.sam" $BAMPATH"/"$BAMNAME".sort.bam" > $BAMPATH"/"$BAMNAME".head.bam"

java $MEM -Djava.io.tmpdir=/tmp -jar ~/bin/picard/ReorderSam.jar I=$BAMPATH"/"$BAMNAME".head.bam" O=$BAMPATH"/"$BAMNAME".lexy.bam" REFERENCE= ~/bin/hg19/ucsc.hg19.fasta CREATE_INDEX=true
