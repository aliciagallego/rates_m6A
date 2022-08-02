#!/usr/bin/env Rscript

####################################
## Retain the first 20 Kb from TSS # 
####################################

# This script tranforms gene coordinates to keep TSS and TSS+20kb
# This gene list is generated to just retain those genes that show transcription from the TSS

# -------
# Paths |
# -------
refseq_path <- "/path/RefSeq_genes.bed"
output <- "/path/output"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=FALSE,
                     col.names = c("Chr","Start","End","Gene_name","NA1","Strand"))

# ----------------------------------------
# Transform start - 2Kb and start + 50Kb |
# ----------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,6] == '+') {
    refseq[i,7]=(refseq[i,2])
    refseq[i,8]=(refseq[i,2]+20000)
  }
  else if(refseq[i,6] == '-'){
    refseq[i,7]=(refseq[i,3]-20000)
    refseq[i,8]=(refseq[i,3])
  }
  else {
    refseq[i,7]=NA
    refseq[i,8]=NA
  }
  count = count+1
}
print(count)

# Replace col 2 by col 7 and col 3 by col 8
count <- 0
for(i in 1:nrow(refseq)) {
  refseq[i,2]=refseq[i,7]
  refseq[i,3]=refseq[i,8]
  count = count+1
}
print(count)

# Drop columns 7 and 8
refseq <- refseq[-c(7,8)]

# -----------------------------------------
# Remove transcripts with negative coords |
# -----------------------------------------
# This step removes mitochondrial and rare sequences
refseq <- refseq[!refseq$Start < 0,] 

# -----------
# Save data |
# -----------
write.table(refseq, 
            file = paste0(output,"RefSeq_TSS_20Kb.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)
