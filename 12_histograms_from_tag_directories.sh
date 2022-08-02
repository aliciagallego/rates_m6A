#!/bin/bash

# This script generates histograms based on the Tag Directories and a list of genes in a BED file (transcriptome)
# The transcriptome contains gene coordinates TSS-2kb and TSS+50kb
# annotatePeaks.pl is a function from HOMER

TAG="/path/TagDirectory/"
TRANSCRIPTOME="/path/RefSeq_coding_TSS_2Kb50Kb.bed
OUT="/path/Histogram"      

annotatePeaks.pl $TRANSCRIPTOME none -size 52000 -hist 20 -ghist -d $TAG > $OUT/output.txt"     	 	
