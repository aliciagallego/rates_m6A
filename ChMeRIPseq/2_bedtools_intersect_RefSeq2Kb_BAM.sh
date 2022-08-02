#!/bin/bash

# This script intersects (bedtools function) the aligned BAM sequences and the BED transcriptome previously modified to include TSS-2kb and TTS+2kb coordinates

AFILE=/path/RefSeq_LongList_2Kb.bed
INPUTS=/path/Alignments/Inputs/*001
IPS=/path/Alignments/IPs/*001

for DIR in $INPUTS
do
	NAME=$(basename $DIR)
	FILE=$DIR/accepted_hits.bam	
	bedtools intersect -a $AFILE -b $FILE -c -s > /path/$NAME.bed
done

for DIR in $IPS
do
	NAME=$(basename $DIR)
	FILE=$DIR/accepted_hits.bam	
	bedtools intersect -a $AFILE -b $FILE -c -s > /path/$NAME.bed	
done
