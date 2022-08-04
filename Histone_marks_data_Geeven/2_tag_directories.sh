#!/bin/bash

# Tag directory is a directory that contains several files describing the sequencing data by chromosome. 
# The input file for making a tag directory are the .sorted.bam files (generated from *.sam aligments)
# makeTagDirectory is a function from HOMER

for HIST in H3K4me1 H3K4me3 H3K9me3 H3K27me3 
do 
	for REP in WT TKO
	do
		for GROUP in slow med fast
		do	
		BAM="/path/"$HIST"/"$HIST"_"$REP"_Pull_sorted.bam"
		OUTPUTDIR="/path/TagDirectories/"$HIST"/"$REP"_"$GROUP"_TagDirectory/" 
		makeTagDirectory $OUTPUTDIR $BAM -format sam 
done
done
done
wait
