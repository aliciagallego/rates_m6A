#!/usr/bin/env Rscript

###################################################################################################################
## Normalize meRIP data (after Ensemble intersection) and merge with elongation rate data from Maslon et al. 2019 #
###################################################################################################################

# This script normalizes meRIP data (m6A levels) by total read number per replicate and IP/Input ratio
# Mean values are calculate for both replicates in WT
# Correlation test between m6A levels and elongation rates reported by Maslon et al. (2019) is computed in WT

# -----------
# Libraries |
# -----------
libraty("ggplot")
library("ggpubr")

# -------
# Paths |
# -------
m6A_path <- "/path/Maslon2019/Intersect_bedtools/"
rates_Maslon_path <- "/path/Maslon2019/Intersect_bedtools/Maslon_data/Elongation_rates_MaslonWT.txt"
rates_out <- "/path/Maslon2019/"

# ---------------------------------------
# Open TXT data Maslon elongation rates |
# ---------------------------------------
rates_Maslon <- read.table(rates_Maslon_path,h=T,sep="\t",stringsAsFactors=FALSE)

# ---------------------------------------------
# Transform Maslon elongation rates to kb/min |
# ---------------------------------------------
rates_Maslon[,4] = rates_Maslon[,4]/1000
colnames(rates_Maslon) <- c("enst_id","R_squared","p_value","Maslon_Kb_min","tx_length","FPKM_control","boundary_5min",
                        "boundary_15min","gene_short_name")

# -------------------------------------
# Open BED data meRIP-Ensembl aligned |
# -------------------------------------
BED_list = list.files(m6A_path, pattern="*.bed")

for (i in seq_along(BED_list)) {
  filename <- sub("_2Kb.bed*", "", BED_list[i])
  df <- read.table(paste0(m6A_path,BED_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene.stable.ID","Transcript.stable.ID","Strand",
                                 "NA1","NA2","Gene_name","Gene_type","meRIP"))
  df <- subset(df, select=-c(NA1,NA2)) # remove uninformative columns
  df <- df[!duplicated(df$Transcript.stable.ID),] # remove duplicated rows
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

# --------------------
# Merge WT m6A files |
# --------------------
mergedWT <- Reduce(function(x,y) merge(x = x, y = y, by = "Transcript.stable.ID", sort = F),
                   list(WU1_IP[,c("Chr","Start","End","Gene.stable.ID","Strand","Gene_name",
                                  "Gene_type","meRIP_norm","Transcript.stable.ID")],
                        WU2_IP[,c("meRIP_norm","Transcript.stable.ID")]))

mergedWT <- mergedWT[,c("Chr","Start","End","Gene.stable.ID","Transcript.stable.ID","Strand","Gene_name",
                        "Gene_type","meRIP_norm.x","meRIP_norm.y")]

colnames(mergedWT) <- c("Chr","Start","End","Gene.stable.ID","Transcript.stable.ID","Strand","Gene_name",
                        "Gene_type","meRIP_WT1","meRIP_WT2")

mergedWT$meRIP_WT <- (mergedWT$meRIP_WT1+mergedWT$meRIP_WT2)/2

# -------------------------------
# Merge WT m6A and Maslon files |
# -------------------------------
rates_Maslon_m6A_WT<- Reduce(function(x,y) merge(x = x, y = y, by.x = "Transcript.stable.ID",by.y = "enst_id", sort = F),
                             list(mergedWT[,c("Chr", "Start", "End", "Gene.stable.ID","Strand", "Gene_name", "Gene_type",
                                              "meRIP_WT1","meRIP_WT2","meRIP_WT","Transcript.stable.ID")],
                                  rates_Maslon[,c("Maslon_Kb_min","gene_short_name","enst_id")]))

lost <- nrow(rates_Maslon) - nrow(rates_Maslon_m6A_WT)
print(paste(lost,"transcripts lost from Maslon list after merge"))

# ---------------------
# Remove m6A outliers |
# ---------------------
Q <- quantile(rates_Maslon_m6A_WT$meRIP_WT, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(rates_Maslon_m6A_WT$meRIP_WT)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_WT <- subset(rates_Maslon_m6A_WT, rates_Maslon_m6A_WT$meRIP_WT > low & rates_Maslon_m6A_WT$meRIP_WT < up)

lost2<-nrow(rates_Maslon)-nrow(saved_values_WT)
print(paste(lost2,"transcripts lost from Maslon list after merge and meRIP outlier removal"))

# -----------
# Save data |
# -----------
write.table(saved_values_WT, 
            file = paste0(rates_out,"meRIP_Maslon/meRIP_Normalized_Maslon_WT_noOutliers_2Kb.txt"),
            quote = F, sep="\t", col.names = T, row.names = F) 

# -------------------------------------------------------------------------
# m6A and elongation rate data from Maslon et al. 2019 - Correlation test |
# -------------------------------------------------------------------------
png(file = paste0(rates_out, "Plots/Pearson_MaslonElongation_meRIP_WT_2Kb.png"))
ggscatter(saved_values_WT, y = "meRIP_WT", x = "Maslon_Kb_min",
          color="lightblue3",shape = 20,
          size=1.5,
          add.params = list(color = "lightseagreen", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "m6A levels", xlab = "Elongation rate Maslon data (Kb/min)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Pearson correlation m6A levels vs. elongation rates (Maslon data) - WT (746 genes + 2Kb)")
dev.off()
