#!/usr/bin/env Rscript

####################################################################
## Define promoter coordinates (TSS±2Kb) from RefSeq transcriptome #
####################################################################

# This script uses a transcriptome BED list as input and defines the promoter coordinates as TSS±2Kb

# -------
# Paths |
# -------
transcriptome_path <- "/path/RefSeq_genes.bed"
groups_path <- "/path/rate_calculation_tertiles/"

# path for output files
output <- "/path/Geeven/Intersect/"

# --------------------------------------------------------
# Read all files from directories and put them in a list |
# --------------------------------------------------------
groups_list = list.files(groups_path, pattern="*txt")

for (i in seq_along(groups_list)) {
  filename <- sub("_without05_", "_", groups_list[i])
  filename <- sub(".txt", "", filename)
  df <- read.table(paste0(groups_path,groups_list[i]), header = T, sep="\t", stringsAsFactors=FALSE)
  assign(filename, df)
}

# -----------
# Open data |
# -----------
transcriptome <- read.table(transcriptome_path,h=F,sep="\t",stringsAsFactors=FALSE,
                     col.names = c("Chr","Start","End","Gene_name","NA1","Strand"))

# -------------------------------------
# Transform start -4Kb and start +4Kb |
# -------------------------------------
count <- 0
for(i in 1:nrow(transcriptome)) {
  if(transcriptome[i,6] == '+') {
    transcriptome[i,7]=(transcriptome[i,2]-2000)
    transcriptome[i,8]=(transcriptome[i,2]+2000)
  }
  else if(transcriptome[i,6] == '-'){
    transcriptome[i,7]=(transcriptome[i,3]-2000)
    transcriptome[i,8]=(transcriptome[i,3]+2000)
  }
  else {
    transcriptome[i,7]=NA
    transcriptome[i,8]=NA
  }
  count = count+1
}
print(count)

# Replace col 2 by col 7 and col 3 by col 8
count <- 0
for(i in 1:nrow(transcriptome)) {
  transcriptome[i,2]=transcriptome[i,7]
  transcriptome[i,3]=transcriptome[i,8]
  count = count+1
}
print(count)

# Drop columns 7 and 8
transcriptome <- transcriptome[-c(7,8)]
head(transcriptome)
tail(transcriptome)

# ------------------------------
# Transcriptome by rate groups |
# ------------------------------
merged_non05 <- merged_noNan[which(merged_noNan$WT_Pull > 0.5 & merged_noNan$TK0_Pull > 0.5),]
refseq$enst_id <- biomart$Transcript.stable.ID[pmatch(refseq$RefseqID, biomart$RefSeq.mRNA.ID, duplicates.ok = FALSE)]
overlap <- refseq[refseq$RefseqID %in% biomart[biomart$Transcript.stable.ID %in% maslon$enst_id, ]$RefSeq.mRNA.ID, ]

transcriptome_WTslow <- transcriptome[transcriptome$Gene_name %in% Slow_WT$Gene_name, ]
transcriptome_WTmed <-  transcriptome[transcriptome$Gene_name %in% Medium_WT$Gene_name, ]
transcriptome_WTfast <- transcriptome[transcriptome$Gene_name %in% Fast_WT$Gene_name, ]
transcriptome_TKOslow <-transcriptome[transcriptome$Gene_name %in% Slow_TKO$Gene_name, ]
transcriptome_TKOmed <- transcriptome[transcriptome$Gene_name %in% Medium_TKO$Gene_name, ]
transcriptome_TKOfast <-transcriptome[transcriptome$Gene_name %in% Fast_TKO$Gene_name, ]

# -----------
# Save data |
# -----------
write.table(transcriptome_WTslow, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TSS_4kb_WT_slow.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_WTmed, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TSS_4kb_WT_med.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_WTfast, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TSS_4kb_WT_fast.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_TKOslow, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TSS_4kb_TKO_slow.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_TKOmed, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TSS_4kb_TKO_med.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_TKOfast, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TSS_4kb_TKO_fast.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)
