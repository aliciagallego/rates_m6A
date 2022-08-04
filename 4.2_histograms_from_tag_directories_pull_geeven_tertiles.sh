#!/bin/bash

for HIST in H3K4me1 H3K4me3 H3K9me3 H3K27me3 

do 
	for REP in WT TKO
	do
		for GROUP in slow med fast
		do
		TAG="/media/cc/B/Alicia/Geeven/Geeven_output/2.2_TagDirectories_tertiles/"$HIST"/"$REP"_"$GROUP"_TagDirectory/"
		OUT="/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/"$HIST"/"$REP"_"$GROUP
		
		annotatePeaks.pl "/media/cc/B/Alicia/Geeven/Geeven_output/1_RefSeq_transcriptome/RefSeq_20621_Transcriptome_TSS_4kb_"$REP"_"$GROUP".bed" none -size 4000 -hist 20 -ghist -d $TAG > $OUT/$HIST"_"$REP"_"$GROUP"output.txt"     	

done
done
done
wait
