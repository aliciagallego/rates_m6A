#!/bin/bash

# This script intersects (bedtools function) the aligned BAM sequences from Xu et al. (2021) and the BED transcriptome 
# previously modified to include TSS±2Kb or TTS±2Kb coordinates

BAMS=/path/Mettl3_Xu/Alignments/Pull/

for FILE in $BAMS/*.bam
do
	for REGION in TSS TTS
	do

	NAME=$(basename $FILE) 
	AFILE="/path/Input_genes_RefSeq_Long_List_2Kb"$REGION".bed"
       	bedtools intersect -a $AFILE -b $FILE -c > /mpath/Mettl3/Intersect/'intersect_'$NAME'_2kb'$REGION.bed

done
done
