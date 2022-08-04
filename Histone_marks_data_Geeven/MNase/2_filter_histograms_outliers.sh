#!/bin/bash

for CEL in WT TKO
do
	for REP in I II
	do
		for GROUP in slow med fast
		do
			INPUT="/path/MNase/Histogram/"$CEL"_"$REP"/"$CEL"_"$REP"_"$GROUP"_output.txt"
			OUT="/path/MNase/Histogram/"$CEL"_"$REP				

			# FRAC value ihas to be defined after prior analyses 
	                # In this case the 0.1% of the most expressed genes will we removed
			FRAC=0.01	

			# get total number of rows (genes) excluding the header 
			MYNUM=$(cat $INPUT | sed '1d' | wc -l)

			# multiply total row number by 0.05 and sum +1
			# sed remove decimals e.g. 4.55 -> 4
			MYNUM=$(echo "($MYNUM * $FRAC) + 1" | bc -l | sed 's/\..*//')

			# Threshold
			# sed "1d" Prints the contents of file excluding the first line to the standard output
			# sum all columns by row
			# sed 's/\..*//':  removes decimals
			# sort -nr: sorts from higher to lower values
			# head -n $MYNUM: selects the first $MYNUM values
			# tail -1: selects the last value from the $MYNUM selectin (it is the value from which the higher values will be removed or threshold)
			THRES=$(cat $INPUT | sed '1d' | awk -v OFS="\t" '{for(i=2;i<=NF;i++) sum+=$i; print sum; sum=0}' | sed 's/\..*//' | sort -nr  | head -n $MYNUM | tail -1)

			printf "up on the threshold\t"$MYNUM"\n" > $OUT/$CEL"_"$REP"_"$GROUP"_summary_"$FRAC".txt"
			printf "threshold\t"$THRES"\n" >> $OUT/$CEL"_"$REP"_"$GROUP"_summary_"$FRAC".txt"

			# prints headers
			cat $INPUT | head -1 > $OUT/$CEL"_"$REP"_"$GROUP"_matrix_"$FRAC".txt"

			# sed "1d": Prints the contents of file excluding the first line to the standard output
			# sum all columns by row
			# selects all rows which summatory column is lower than the previously computed threshold
			cat $INPUT | sed '1d' | awk -v OFS="\t" '{for(i=2;i<=NF;i++) sum+=$i; print $0,sum; sum=0}' | awk -v OFS="\t" -v thres=$THRES '($NF<thres){printf $1"\t"; for (i=2; i<NF; i++) printf $i"\t"; printf "\n"}' >> $OUT/$CEL"_"$REP"_"$GROUP"_matrix_"$FRAC".txt"     	 	

done
done
done
wait
