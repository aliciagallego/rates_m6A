#!/bin/bash

for HIST in H3K4me1 H3K4me3 H3K9me3 H3K27me3 
do 
	for REP in WT TKO
	do
		for GROUP in slow med fast
		do 

			FRAC=0.01	
			INPUT="/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/"$HIST"/"$REP"_"$GROUP"/"$HIST"_"$REP"_"$GROUP"_matrix_"$FRAC".txt"  
			OUT="/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/"$HIST"/"$REP"_"$GROUP

			cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'> $OUT/$HIST"_"$REP"_"$GROUP"_matrix_sum_"$FRAC".txt"

done
done
done


