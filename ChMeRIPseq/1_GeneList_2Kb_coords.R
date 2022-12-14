#!/usr/bin/env Rscript

##########################################################
## Include 2Kb region at TSS (start) and TTS (end) sides #
##########################################################

# This script uses a transcriptome BED list as input and modifies gene coordinates to include at both sides of the gene body

# -------
# Paths |
# -------
refseq_path <- "/path/Genome_files/RefSeq_genes.bed"
output <- "/path/meRIP/Gene_selection_RefSeq/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=FALSE)

# -------------------------------------
# Transform start - 2Kb and end + 2Kb |
# -------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,6] == '+' | refseq[i,6] == '-') {
    refseq[i,2]=(refseq[i,2]-2000)
    refseq[i,3]=(refseq[i,3]+2000)
  }
  else {
    refseq[i,8]=NA
    refseq[i,9]=NA
  }
  count = count+1
}
print(count)

# ----------------------------------------------------------------------------
# Remove transcripts with negative coords (mitochondrial and rare sequences) |
# ----------------------------------------------------------------------------
refseq <- refseq[!refseq$V2 < 0,]

# -----------
# Save data |
# -----------
write.table(refseq, 
            file = paste0(output,"RefSeq_LongList_2Kb.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)
