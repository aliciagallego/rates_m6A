#!/bin/bash

for CEL in WT TKO
do
	for REP in I II
	do
		for GROUP in slow med fast
		do
		
		BED="/media/cc/B/Alicia/Geeven/Geeven_output/1_RefSeq_transcriptome/RefSeq_20621_Transcriptome_TSS_4kb_"$CEL"_"$GROUP".bed"	
		BEDGRAPH="/media/cc/C/Alicia/MNase/bedGraphs/MN_H1_"$CEL"_"$REP"_sorted_depth_wl_trimmed_PE.bedGraph"
	
		OUT="/media/cc/C/Alicia/MNase/histograms/"$CEL"_"$REP

		annotatePeaks.pl $BED none -size 4000 -hist 20 -ghist -bedGraph $BEDGRAPH > $OUT/$CEL"_"$REP"_"$GROUP"_output.txt"     	 	


done
done
done
wait
