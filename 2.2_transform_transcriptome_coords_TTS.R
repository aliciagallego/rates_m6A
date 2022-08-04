#!/usr/bin/env Rscript

#####################
## Include TTS-+4Kb #
#####################

# -------
# Paths |
# -------
# de novo assembly transcriptome from Josemi (cheRNA data) just containing protein coding genes
transcriptome_path <- "/media/cc/A/Alicia/Genome_files/Josemi/RefSeq_genes.bed"
groups_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/2_Rate_calculation_tertiles/"

# path for output files
output <- "/media/cc/B/Alicia/Geeven/Geeven_output/1_RefSeq_transcriptome/"

# --------------------------------------------------------
# Read all files from directories and put them in a list |
# --------------------------------------------------------
groups_list = list.files(groups_path, pattern="*txt")
groups_list

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
head(transcriptome)
nrow(transcriptome)

# -------------------------------------
# Transform start -4Kb and start +4Kb |
# -------------------------------------
count <- 0
for(i in 1:nrow(transcriptome)) {
  if(transcriptome[i,6] == '+') {
    transcriptome[i,7]=(transcriptome[i,3]-2000)
    transcriptome[i,8]=(transcriptome[i,3]+2000)
  }
  else if(transcriptome[i,6] == '-'){
    transcriptome[i,7]=(transcriptome[i,2]-2000)
    transcriptome[i,8]=(transcriptome[i,2]+2000)
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

# -------------
# Rate groups |
# -------------

# -----------
# Save data |
# -----------
write.table(transcriptome, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TSS_4kb.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

# ------------------------------
# Transcriptome by rate groups |
# ------------------------------
nrow(transcriptome)
nrow(Fast_TKO)

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
            file = paste0(output,"RefSeq_20621_Transcriptome_TTS_4kb_WT_slow.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_WTmed, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TTS_4kb_WT_med.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_WTfast, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TTS_4kb_WT_fast.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_TKOslow, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TTS_4kb_TKO_slow.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_TKOmed, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TTS_4kb_TKO_med.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)

write.table(transcriptome_TKOfast, 
            file = paste0(output,"RefSeq_20621_Transcriptome_TTS_4kb_TKO_fast.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)