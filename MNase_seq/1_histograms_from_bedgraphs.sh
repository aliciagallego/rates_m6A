#!/bin/bash

# This script generates histograms using bedGraph files from aligned data as input
# The transcriptome contains gene coordinates TSS±2kb
# annotatePeaks.pl is a function from HOMER

for CEL in WT TKO
do
	for REP in I II
	do
		for GROUP in slow med fast
		do
		
		BED="/path/RefSeq_TSS2kb_"$CEL"_"$GROUP".bed"	
		BEDGRAPH="/path/bedGraphs/"$CEL"_"$REP".bedGraph"
		OUT="/path/MNase/Histogram/"$CEL"_"$REP
		annotatePeaks.pl $BED none -size 4000 -hist 20 -ghist -bedGraph $BEDGRAPH > $OUT/$CEL"_"$REP"_"$GROUP"_output.txt"     	 	


done
done
done
wait
