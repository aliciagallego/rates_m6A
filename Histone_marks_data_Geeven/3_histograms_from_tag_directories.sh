#!/bin/bash

# This script generates histograms based on the Tag Directories and a list of genes in a BED file (transcriptome)
# The transcriptome contains gene coordinates TSS±2kb
# annotatePeaks.pl is a function from HOMER

for HIST in H3K4me1 H3K4me3 H3K9me3 H3K27me3 

do 
	for REP in WT TKO
	do
		for GROUP in slow med fast
		do
		TAG="/path/TagDirectories/"$HIST"/"$REP"_"$GROUP"_TagDirectory/"
		OUT="/path/Histograms/"$HIST"/"$REP"_"$GROUP
		annotatePeaks.pl "path/RefSeq_Transcriptome_TSS2kb_"$REP"_"$GROUP".bed" none -size 4000 -hist 20 -ghist -d $TAG > $OUT/$HIST"_"$REP"_"$GROUP"output.txt"     	

done
done
done
wait
