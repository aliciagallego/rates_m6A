#!/usr/bin/env Rscript

##########################
## Normalize cheRNA data #
##########################

# This script normalizes cheRNA data (chromatin enriched RNA expression levels) by total read number per replicate and transcript size
# Mean values are calculate for both replicates in WT and H1-TKO samples
# Boxplots are computed for both WT and H1-TKO

# -------
# Paths |
# -------
cheRNA_path <- "/path/cheRNA/Intersect_RefSeq/"
output <- "/path/cheRNA/Normalized/"

# ------------
# Open files |
# ------------
cheRNA_list = list.files(cheRNA_path, pattern="*.bed")

for (i in seq_along(cheRNA_list)) {
  filename <- sub("_RefSeq_intersect.bed", "", cheRNA_list[i])
  df <- read.table(paste0(cheRNA_path,cheRNA_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene_name", "NA1", "Strand",
                                 "cheRNA"))
  df <- subset(df, select=-c(NA1)) # remove uninformative columns
  assign(filename, df)
}

# ------------------------------------
# Total reads (data from experiment) |
# ------------------------------------
TKO_I_reads	<-   93717203
TKO_II_reads	<- 93903878
TKO_III_reads	<- 98738497
WT_I_reads	<-   85262714
WT_II_reads	<-   98517689
WT_III_reads	<- 83632131

# --------------------------------------------------
# Normalization by total reads and transcript size |
# --------------------------------------------------

WT_I$cheRNA_normI = ((WT_I$cheRNA/WT_I_reads) / ((WT_I$End-WT_I$Start)/1000)) 
WT_II$cheRNA_normII = ((WT_II$cheRNA/WT_II_reads) / ((WT_II$End-WT_II$Start)/1000)) 
WT_III$cheRNA_normIII = ((WT_III$cheRNA/WT_III_reads) / ((WT_III$End-WT_III$Start)/1000)) 

TKO_I$cheRNA_normI = ((TKO_I$cheRNA/TKO_I_reads) / ((TKO_I$End-TKO_I$Start)/1000)) 
TKO_II$cheRNA_normII = ((TKO_II$cheRNA/TKO_II_reads) / ((TKO_II$End-TKO_II$Start)/1000)) 
TKO_III$cheRNA_normIII = ((TKO_III$cheRNA/TKO_III_reads) / ((TKO_III$End-TKO_III$Start)/1000)) 

# ------------------------------
# Get means between replicates |
# ------------------------------
WT_table <- subset(WT_I, select = -c(cheRNA))
TKO_table <- subset(TKO_I, select = -c(cheRNA))

WT_table$cheRNA_normII  = WT_II$cheRNA_normII
WT_table$cheRNA_normIII = WT_III$cheRNA_normIII
TKO_table$cheRNA_normII = TKO_II$cheRNA_normII
TKO_table$cheRNA_normIII= TKO_III$cheRNA_normIII

WT_table$cheRNA_norm_mean  = (WT_table$cheRNA_normI + WT_table$cheRNA_normII + WT_table$cheRNA_normIII)/3
TKO_table$cheRNA_norm_mean = (TKO_table$cheRNA_normI + TKO_table$cheRNA_normII + TKO_table$cheRNA_normIII)/3

# ------------
# Merge data |
# ------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                  list(WT_table,
                       TKO_table[,c("cheRNA_normI","cheRNA_normII","cheRNA_normIII","cheRNA_norm_mean","Gene_name")]))

colnames(merged) <- c("Gene_name", "Chr","Start","End","Strand",
                            "cheRNA_WT1","cheRNA_WT2","cheRNA_WT3","cheRNA_WT_mean",
                            "cheRNA_TKO1","cheRNA_TKO2","cheRNA_TKO3","cheRNA_TKO_mean")

# -----------
# Save data |
# -----------
write.table(merged, 
            file = paste0(output,"cheRNA_RefSeq_LongList_Normalized.bed"),
            quote = F, sep="\t", col.names = T, row.names = F)     

# ----------
# Boxplots |
# ----------
wilcox <- wilcox.test(merged$cheRNA_WT_mean, merged$cheRNA_TKO_mean)
png(file = paste0(output_plots, "BoxPlots_cheRNAs_means.png"))
boxplot(merged$cheRNA_WT_mean, merged$cheRNA_TKO_mean,
        col = c(4,2),
        at = c(1,2),
        names = c("WT", "TKO"),
        main = "cheRNA (means) (20621 genes)",
        ylab="cheRNA reads",
        outline = F,
        boxwex = 0.3)
dev.off()	
