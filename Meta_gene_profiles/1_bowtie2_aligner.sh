#!/bin/bash

READS1="/path/DRB-4sU/R1_001.fastq.gz"
READS2="/path/DRB-4sU/R2_001.fastq.gz"
GENOME="path/Genome_files/mm10/Bowtie_index"
OUTPUT="/path/alignment.sam"

bowtie2 -x $GENOME -1 $READS1 -2 $READS2 -S $SALIDA -N 1
