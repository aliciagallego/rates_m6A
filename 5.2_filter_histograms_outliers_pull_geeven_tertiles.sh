#!/bin/bash

for HIST in H3K4me1 H3K4me3 H3K9me3 H3K27me3 
do 
	for REP in WT TKO
	do

		for GROUP in slow med fast
		do
			INPUT="/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/"$HIST"/"$REP"_"$GROUP"/"$HIST"_"$REP"_"$GROUP"output.txt" 
			OUT="/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/"$HIST"/"$REP"_"$GROUP
			FRAC=0.01	

			# obtener el numero total de filas (genes) sin contar la fila del encabezado
			MYNUM=$(cat $INPUT | sed '1d' | wc -l)

			# multiplicar el numero total de filas por $FRAC y sumar 1
			# bc -l indica usar formulas matemáticas
			# sed quita los decimales 4.55 -> 4
			MYNUM=$(echo "($MYNUM * $FRAC) + 1" | bc -l | sed 's/\..*//')

			# Threshold
			# sed "1d" Prints the contents of file excluding the first line to the standard output
			# sumar todas las columnas por fila
			# sed 's/\..*//':  quita los decimales
			# sort -nr: ordena de mayor a menor
			# head -n $MYNUM: selecciona los $MYNUM primeros
			# tail -1: selecciona el último de los valores $MYNUM (el número a partir del cual se van a eliminar los valores altos o threshold)

			THRES=$(cat $INPUT | sed '1d' | awk -v OFS="\t" '{for(i=2;i<=NF;i++) sum+=$i; print sum; sum=0}' | sed 's/\..*//' | sort -nr  | head -n $MYNUM | tail -1)

			printf "up on the threshold\t"$MYNUM"\n" > $OUT/$HIST"_"$REP"_"$GROUP"_summary_"$FRAC".txt"
			printf "threshold\t"$THRES"\n" >> $OUT/$HIST"_"$REP"_"$GROUP"_summary_"$FRAC".txt"

			# imprimir los encabezados
			cat $INPUT | head -1 > $OUT/$HIST"_"$REP"_"$GROUP"_matrix_"$FRAC".txt"

			# sed "1d": Prints the contents of file excluding the first line to the standard output
			# sumar todas las columnas por fila
			# seleccionar aquellas filas cuya columna sumatoria es menor al threhold anteriormente calculado
			cat $INPUT | sed '1d' | awk -v OFS="\t" '{for(i=2;i<=NF;i++) sum+=$i; print $0,sum; sum=0}' | awk -v OFS="\t" -v thres=$THRES '($NF<thres){printf $1"\t"; for (i=2; i<NF; i++) printf $i"\t"; printf "\n"}' >> $OUT/$HIST"_"$REP"_"$GROUP"_matrix_"$FRAC".txt"     	 	

done
done
done
wait
