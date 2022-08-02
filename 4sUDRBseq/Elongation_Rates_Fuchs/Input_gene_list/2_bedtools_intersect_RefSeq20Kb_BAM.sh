#!/bin/bash

# This script intersects (bedtools function) the aligned BAM sequences and the BED transcriptome previously modified to just include TSS and TSS+20kb coordinates
# The intersection is just performed over the 5 min 4sU-DRB samples

for NUM in 1 2 3 4
do
	BAM="/path/4sUDRB/MG9-"$NUM"_filtered.bam"
	BED="/path/RefSeq_LongList_TSS_20Kb.bed"
	SALIDA="/path/Intersect_RefSeqBED20Kb_BAM/"$NUM".bed"

	bedtools intersect -a $BED -b $BAM -c -s > $SALIDA
done
