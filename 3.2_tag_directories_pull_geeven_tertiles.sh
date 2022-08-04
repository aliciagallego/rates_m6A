#!/bin/bash

# Tag directory: is a directory that contains several files describing your data. In this case the input file for making a tag directory are the .sorted.bam files from the Pulls

for HIST in H3K4me1 H3K4me3 H3K9me3 H3K27me3 
do 
	for REP in WT TKO
	do
		for GROUP in slow med fast
		do
	
		BAM="/media/cc/B/Alicia/Geeven/Geeven_data/"$HIST"/"$HIST"_"$REP"_Pull_sorted.bam"
		OUTPUTDIR="/media/cc/B/Alicia/Geeven/Geeven_output/2.2_TagDirectories_tertiles/"$HIST"/"$REP"_"$GROUP"_TagDirectory/" 
 
		# since we are working with ChIP seq data (histone modifications) we remove -flip -sspe parameters
		makeTagDirectory $OUTPUTDIR $BAM -format sam 
done
done
done
wait
