#!/bin/bash

# This script intersects (bedtools function) the aligned BAM sequences and the BED transcriptome previously modified to include TSS±2Kb coordinates

AFILE=/path/IRefSeq_LongList_2KbTSS.bed
BAMS=/path/H1_Cao/Alignments/

for FILE in $BAMS/*.bam

do
	NAME=$(basename $FILE) 
  bedtools intersect -a $AFILE -b $FILE -c > path/$NAME'.bed

done
