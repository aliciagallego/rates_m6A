# Tag directory is a directory that contains several files describing the sequencing data by chromosome. 
# The input file for making a tag directory are the .sorted.bam files (generated from *.sam aligments)
# makeTagDirectory is a function from HOMER

INPUT="/path/*alinment_sorted.bam"
OUTPUTDIR="/path/TagDirectory/"        	

makeTagDirectory $OUTPUTDIR -flip -sspe $INPUT -format sam 
