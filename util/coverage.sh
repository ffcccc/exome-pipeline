#!/bin/bash
BAM=$1
BED=$2
# e.g ~/bin/hg19/target.intervals.bed

bedtools coverage -hist -abam ${BAM} -b ${BED} | grep ^all > ${BAM}.hist.coverage
