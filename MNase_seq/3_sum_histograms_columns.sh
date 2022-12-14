#!/bin/bash

# This script computes total reads per histogram window across all genes

for CEL in WT TKO
do
	for REP in I II
	do
		for GROUP in slow med fast
		do

			FRAC=0.01	
			INPUT="/path/MNase/Hstogram/"$CEL"_"$REP"/"$CEL"_"$REP"_"$GROUP"_matrix_"$FRAC".txt"
			OUT="/path/MNase/Hstogram/"$CEL"_"$REP

			cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'> $OUT/$CEL"_"$REP"_"$GROUP"_matrix_sum_"$FRAC".txt"

done
done
done
