#!/bin/bash

# This script intersects (bedtools function) the aligned BAM sequences and the BED transcriptome (Ensembl)
# previously modified to include TSS-2kb and TTS+2kb coordinates

AFILE=/path/Maslon2019/Genome_data/Ensembl_mm10_transcripts_2Kb.bed
INPUTS=/path/meRIP/Alignments/Inputs/*001
IPS=/path/meRIP/Alignments/IPs/*001

for DIR in $INPUTS
do
	NAME=$(basename $DIR)
	FILE=$DIR/accepted_hits.bam	
	bedtools intersect -a $AFILE -b $FILE -c -s > /path/Maslon2019/Intersect_bedtools/$NAME.bed
done

for DIR in $IPS
do
	NAME=$(basename $DIR)
	FILE=$DIR/accepted_hits.bam	
	bedtools intersect -a $AFILE -b $FILE -c -s > /path/Maslon2019/Intersect_bedtools/$NAME.bed
done
