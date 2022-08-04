#!/bin/bash

# This script generates histograms based on the Tag Directories and a list of genes in a BED file (transcriptome)
# The transcriptome contains gene coordinates TSS-2kb and TSS+50kb
# annotatePeaks.pl is a function from HOMER

TRANSCRIPTOME="/path/RefSeq_coding_TSS_2Kb50Kb.bed"

for NUM in WT0 WT5 WT15 WT45 TKO0 TKO5 TKO15 TKO45
do

	TAG="/path/"$NUM"_TagDirectory/"
	OUT="/path/Histograms/"$NUM"_Histogram"
	annotatePeaks.pl $TRANSCRIPTOME none -size 52000 -hist 20 -ghist -d $TAG > $OUT/"$NUM"_output.txt"     	 	


done
wait
