#!/usr/bin/env Rscript

#########################
## Normalize meRIP data #
#########################

# This script normalizes meRIP data (m6A levels) by total read number per replicate and IP/Input ratio
# Mean values are calculate for both replicates in WT and H1-TKO samples
# Boxplots are computed for both WT and H1-TKO

# -------
# Paths |
# -------
m6A_path <- "/path/meRIP/Intersect_RefSeq_BAM/"
output <- "/path/meRIP/Normalized/"

# ---------------
# Open BED data |
# ---------------
BED_list = list.files(m6A_path, pattern="*.bed")

for (i in seq_along(BED_list)) {
  filename <- sub("_2Kb.bed*", "", BED_list[i])
  df <- read.table(paste0(m6A_path,BED_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene_name","NA1","Strand","meRIP"))
  df <- subset(df, select=-c(NA1)) # remove uninformative columns
  df <- df[!duplicated(df$Gene_name),] # remove duplicated rows
  assign(filename, df)
}

# ------------------------------------
# Total reads (data from experiment) |
# ------------------------------------
WU1_IP_reads	<- 52603601
WU2_IP_reads	<- 56468821
TU1_IP_reads	<- 52209521
TU2_IP_reads	<- 53790404

WU1_In_reads	<- 35203047
WU2_In_reads	<- 36358279
TU1_In_reads	<- 34963985
TU2_In_reads	<- 34564725

# ------------------------------------
# Normalize by total reads and IP/In |
# ------------------------------------
WU1_IP$meRIP_norm = (WU1_IP$meRIP/WU1_IP_reads) / (WU1_In$meRIP/WU1_In_reads)
WU2_IP$meRIP_norm = (WU2_IP$meRIP/WU2_IP_reads) / (WU2_In$meRIP/WU2_In_reads)
TU1_IP$meRIP_norm = (TU1_IP$meRIP/TU1_IP_reads) / (TU1_In$meRIP/TU1_In_reads)
TU2_IP$meRIP_norm = (TU2_IP$meRIP/TU2_IP_reads) / (TU2_In$meRIP/TU2_In_reads)

# -----------------------------
# Remove NA, Inf and 0 values |
# -----------------------------
WU1_IP <- subset(WU1_IP, meRIP_norm!="NaN" & meRIP_norm!=0 & meRIP_norm!="Inf")
WU2_IP <- subset(WU2_IP, meRIP_norm!="NaN" & meRIP_norm!=0 & meRIP_norm!="Inf")
TU1_IP <- subset(TU1_IP, meRIP_norm!="NaN" & meRIP_norm!=0 & meRIP_norm!="Inf")
TU2_IP <- subset(TU2_IP, meRIP_norm!="NaN" & meRIP_norm!=0 & meRIP_norm!="Inf")

# -----------------
# Merge m6A files |
# -----------------

# WT
# --
mergedWT <- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                   list(WU1_IP[,c("Chr","Start","End","Gene_name","Strand","meRIP_norm")],
                        WU2_IP[,c("meRIP_norm","Gene_name")]))

colnames(mergedWT) <- c("Gene_name","Chr","Start","End","Strand","meRIP_WT1","meRIP_WT2")
mergedWT$meRIP_WT <- (mergedWT$meRIP_WT1+mergedWT$meRIP_WT2)/2


# TKO
# ---
mergedTKO <- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                    list(TU1_IP[,c("Chr","Start","End","Gene_name","Strand","meRIP_norm")],
                         TU2_IP[,c("meRIP_norm","Gene_name")]))

colnames(mergedTKO) <- c("Gene_name","Chr","Start","End","Strand","meRIP_TKO1","meRIP_TKO2")
mergedTKO$meRIP_TKO <- (mergedTKO$meRIP_TKO1+mergedTKO$meRIP_TKO2)/2


# -----------
# Save data |
# -----------
write.table(mergedWT, 
            file = paste0(output,"Normalized_data/meRIP_RefSeqLongList_Normalized_WT_2Kb.txt"),
            quote = F, sep="\t", col.names = T, row.names = F) 

write.table(mergedTKO, 
            file = paste0(output,"Normalized_data/meRIP_RefSeqLongList_Normalized_TKO_2Kb.txt"),
            quote = F, sep="\t", col.names = T, row.names = F) 

# ----------------------------
# Mann-Whitney-Wilcoxon test |
# ----------------------------
sink(paste0(output,"4_Plots/meRIP_RefSeqLongList_WilcoxTest.txt"))
wilcox.test(mergedWT$meRIP_WT, mergedTKO$meRIP_TKO)
sink()

# ----------
# Boxplots |
# ----------

# Means
png(file = paste0(output, "4_Plots/meRIP_All_RefSeqLongList_2Kb"))
boxplot(mergedWT$meRIP_WT, mergedTKO$meRIP_TKO,
        col = c(4,2),
        at = c(1,2),
        names = c("WT", "TKO"),
        main = "m6A - means",
        ylab="m6A levels",
        outline = F,
        boxwex = 0.3)
dev.off()	

# Replicates
png(file = paste0(output, "4_Plots/meRIP_replicates_RefSeqLongList_2Kb"))
boxplot(mergedWT$meRIP_WT1, mergedWT$meRIP_WT2, mergedTKO$meRIP_TKO1,mergedTKO$meRIP_TKO2,
        col = c(4,4,2,2),
        at = c(1,2,4,5),
        names = c("WT1","WT2","TKO1", "TKO2"),
        main = "m6A - replicates",
        ylab="m6A levels",
        outline = F,
        boxwex = 0.5)
dev.off()
