#!/bin/bash

for NUM in 2 3 4 5 7 8 9 10
do

	TAG="/media/cc/B/Josemi/TTseq_LowCoverage/TTseq_LowCoverage_output/Pull_Metaplots/MG9-"$NUM"_TagDirectory/"

	OUT="/media/cc/B/Josemi/TTseq_LowCoverage/TTseq_LowCoverage_output/Pull_Metaplots/Histograms/MG9-"$NUM"_Histogram"

	annotatePeaks.pl /media/cc/B/Josemi/TTseq_Feb2022/TTseq_data/RefSeq_LongList_TSS_2Kb50Kb_OK.bed none -size 52000 -hist 20 -ghist -d $TAG > $OUT/"MG9-"$NUM"_output.txt"     	 	


done
wait
