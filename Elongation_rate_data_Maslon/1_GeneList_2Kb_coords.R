#!/usr/bin/env Rscript

#############################################################
## Include 2Kb region at TSS (start) and TTS (end) sides #
#############################################################

# This script uses a transcriptome BED list (Ensembl) as input and modifies gene coordinates to include at both sides of the gene body

# -------
# Paths |
# -------
ensembl_path <- "/path/Maslon2019/Genome_data/Ensembl_mm10_transcripts.bed"
rates_out    <- "/path/Maslon2019/Genome_data/"

# -----------
# Open data |
# -----------
ensembl <- read.table(ensembl_path,h=F,sep="\t",stringsAsFactors=FALSE)

# ---------------------------------------------------
# Transform coords: start = start-2Kb end = end+2Kb |
# ---------------------------------------------------
count <- 0
for(i in 1:nrow(ensembl)) {
  if(ensembl[i,6] == '+' | ensembl[i,6] == '-') {
    ensembl[i,2]=(ensembl[i,2]-2000)
    ensembl[i,3]=(ensembl[i,3]+2000)
  }
  else {
    ensembl[i,8]=NA
    ensembl[i,9]=NA
  }
  count = count+1
}
print(count)

# ----------------------------------------------------------------------------
# Remove transcripts with negative coords (mitochondrial and rare sequences) |
# ----------------------------------------------------------------------------
ensembl2 <- ensembl[!ensembl$V2 < 0,]

# -----------
# Save data |
# -----------
write.table(ensembl2, 
            file = paste0(rates_out,"Ensembl_mm10_transcripts_2Kb.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)
