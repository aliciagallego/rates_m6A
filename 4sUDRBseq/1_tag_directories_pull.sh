#!/bin/bash

# Tag directory: is a directory that contains several files describing your data. In this case the input file for making a tag directory are the .sorted.bam files (obtained from *.sam)

# makeTagDirectory is a function from HOMER

BAM=/media/cc/A/Josemi/NGS/TTseq/Alignments/Pull/MG9-10_S10

for DIR in $BAM 

do
	NAME=$(basename $DIR | sed 's/_S.*//')
	FILE=$DIR/accepted_hits.bam
	OUTPUTDIR="/media/cc/B/Josemi/TTseq_LowCoverage/TTseq_LowCoverage_output/Pull_Metaplots/"$NAME"_TagDirectory/"  
	
	makeTagDirectory $OUTPUTDIR -flip -sspe $FILE -format sam 

done
wait


	
