#!/bin/bash

for NUM in 2 3 4 5 7 8 9 10

do
	FRAC=0.01	
	INPUT="/media/cc/B/Josemi/TTseq_LowCoverage/TTseq_LowCoverage_output/Pull_Metaplots/Histograms/MG9-"$NUM"_Histogram/MG9-"$NUM"_matrix_"$FRAC".txt"
	OUT="/media/cc/B/Josemi/TTseq_LowCoverage/TTseq_LowCoverage_output/Pull_Metaplots/Histograms/MG9-"$NUM"_Histogram"
	OUT2="/media/cc/B/Josemi/TTseq_LowCoverage/TTseq_LowCoverage_output/Pull_Metaplots/Histograms"

	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'> $OUT/"MG9-"$NUM"_matrix_sum_"$FRAC".txt"

	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'>> $OUT2/"sum_all_matrix_sum_"$FRAC".txt"

done
wait

