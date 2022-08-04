#!/bin/bash

# This script computes total reads per histogram window across all genes

for NUM in WT0 WT5 WT15 WT45 TKO0 TKO5 TKO15 TKO45

do
	FRAC=0.01	
	INPUT="path/Histograms/"$NUM"_Histogram/MG9-"$NUM"_matrix_"$FRAC".txt"
	OUT="/path/Histograms/"$NUM"_Histogram"
	OUT2="/path/Histograms"

	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'> $OUT/"$NUM"_matrix_sum_"$FRAC".txt"
	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'>> $OUT2/"sum_all_matrix_sum_"$FRAC".txt"
done
wait
