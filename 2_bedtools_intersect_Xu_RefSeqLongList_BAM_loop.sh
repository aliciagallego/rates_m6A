#!/bin/bash

BAMS=/media/cc/A/Josemi/Mettl3_ChIP_Xu/Alignments/Pull/

for FILE in $BAMS/*.bam
do
	for REGION in TSS TTS
	do

	NAME=$(basename $FILE | sed 's/.bam//') 
	# bed files containing TSS+-2Kb and TTS +-2Kb also used for Cao et al. on H1c Hid analyses
	AFILE="/media/cc/B/Alicia/H1_Cao/H1_Cao_output/1_bedtools_intersect/Input_genes_RefSeq_Long_List_2Kb"$REGION".bed"
       	bedtools intersect -a $AFILE -b $FILE -c > /media/cc/C/Alicia/Mettl3/Mettl3_output/1_bedtools_intersect/'intersect_'$NAME'_2kb'$REGION.bed

done
done

