#!/usr/bin/env Rscript

##################################
## RNAPII pausing index analysis #
##################################

# This script computes RNAPII pausing index as the ratio between RNAPII occupancies in the promoter and the gene body
# Boxplots are computed for both WT and H1-TKO

# -------
# Paths |
# -------
RNAPII_GB_path <- "/path/RNAPII/Intersect_bedtools/GeneBody/"
RNAPII_PR_path <- "/path/RNAPII/Intersect_bedtools/Promoters/"

output <- "/path/RNAPII/Pausing_index/"
output_plots <- "/path/RNAPII/Pausing_index/Plots/"

# -----------
# Open data |
# -----------
RNAPII_GB <- read.table(RNAPII_GB_path ,h=T, sep="\t", stringsAsFactors=F)
RNAPII_PR <- read.table(RNAPII_PR_path ,h=T, sep="\t", stringsAsFactors=F)

# -------------------
# Merge RNAPII data |
# -------------------
mergedRNAPII_GB_PR<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                         list(RNAPII_GB[,c("Gene_name","Chr","Start","End","Strand","RNAPII_WTmean_GB","RNAPII_TKOmean_GB")],
                              RNAPII_PR[,c("Gene_name","Start","End","RNAPII_WTmean_PR","RNAPII_TKOmean_PR")]))

colnames(mergedRNAPII_GB_PR) <- c("Gene_name","Chr","Start_GB","End_GB","Strand","RNAPII_WT_GB","RNAPII_TKO_GB",
                               "Start_PR","End_PR","RNAPII_WT_PR","RNAPII_TKO_PR")

# -------------------------
# Calculate pausing index |
# -------------------------
mergedRNAPII_GB_PR$PausIndex_WT = mergedRNAPII_GB_PR$RNAPII_WT_PR/mergedRNAPII_GB_PR$RNAPII_WT_GB
mergedRNAPII_GB_PR$PausIndex_TKO = mergedRNAPII_GB_PR$RNAPII_TKO_PR/mergedRNAPII_GB_PR$RNAPII_TKO_GB

# -----------------------------
# Remove NA, Inf and 0 values |
# -----------------------------
mergedRNAPII_GB_PR <- subset(mergedRNAPII_GB_PR, PausIndex_WT!="NaN" & PausIndex_WT!=0 & PausIndex_WT!="Inf")
mergedRNAPII_GB_PR <- subset(mergedRNAPII_GB_PR, PausIndex_TKO!="NaN" & PausIndex_TKO!=0 & PausIndex_TKO!="Inf")

# -----------
# Save data |
# -----------
write.table(mergedRNAPII_GB_PR, 
            file = paste0(output,"RNAPII_RefSeq_LongList_PausingIndex_500bp.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)     

# -------------------------------
# Boxplots RNAPII pausing index |
# -------------------------------
meanwT <- round(mean(mergedRNAPII_GB_PR$PausIndex_WT),digits=2)
meanTKO <- round(mean(mergedRNAPII_GB_PR$PausIndex_TKO), digits=2)
means <- c(meanwT,meanTKO)

png(file = paste0(output_plots, "BoxPlots_RNAPII_RefSeq_LongList_PausingIndex_500bp.png"))
boxplot(mergedRNAPII_GB_PR$PausIndex_WT, mergedRNAPII_GB_PR$PausIndex_TKO,
        col = c(4,2),
        at = c(1,2),
        names = c("WT", "TKO"),
        main = "RNApolII pausing index (means) (20217 genes)",
        ylab="Promoter pausing index",
        outline = F,
        boxwex = 0.3)
text(x=c(1:2),y=4.2,
     paste("mean = ",means,sep=""))
dev.off()	

# -------------------------------------------------
# Mann-Whitney-Wilcoxon test RNAPII pausing index |
# -------------------------------------------------
sink(paste0(output_plots, "RNApolII_WilcoxTest_RefSeq_LongList_PausingIndex_500bp.txt"))
wilcox.test(mergedRNAPII_GB_PR$PausIndex_WT, mergedRNAPII_GB_PR$PausIndex_TKO)
sink()
