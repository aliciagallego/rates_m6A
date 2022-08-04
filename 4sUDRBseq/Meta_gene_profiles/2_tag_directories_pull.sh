# Tag directory is a directory that contains several files describing the sequencing data by chromosome. 
# The input file for making a tag directory are the .sorted.bam files (generated from *.sam aligments)
# makeTagDirectory is a function from HOMER

BAM=/path/Alignments

for DIR in $BAM 
do
	NAME=$(basename $DIR)
	FILE=$DIR/accepted_hits.bam
	OUTPUTDIR="/path/"$NAME"_TagDirectory/"  
	makeTagDirectory $OUTPUTDIR -flip -sspe $FILE -format sam 

done
wait
