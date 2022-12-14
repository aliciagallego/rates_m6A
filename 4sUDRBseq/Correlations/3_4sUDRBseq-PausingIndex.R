#!/usr/bin/env Rscript

#########################################################################
## Correlations elongation rates (4sUDRBseq or TTseq) and Pausing Index #
#########################################################################

# -----------
# Libraries |
# -----------
library("ggplot2")
library("ggpubr")
library("reshape2")

# -------
# Paths |
# -------
Pausing_path <- "/path/RNAPII/Pausing_index/RNAPII_RefSeq_LongList_PausingIndex_500bp.txt"
TTseq_path <- "/path/4sUDRB/Elongation_rate/Rate_calculation/Elongation_rate_20Kb_Pull_processed_without05.txt"

output <- "/path/4sUDRB/Elongation_rate/Correlations/"

# -----------
# Open data |
# -----------
Pausing <- read.table(Pausing_path,h=T,sep="\t",stringsAsFactors=F)
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=F,
                    col.names = c("Gene_name","Chr", "Start","End","Exon_number","Strand","TTseq_WT_pull","TTseq_TKO_pull","Size"))
# -------------
# Merge files |
# -------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                list(Pausing[,c("Chr","Start_GB","End_GB","Strand","PausIndex_WT","PausIndex_TKO","Gene_name")],
                     TTseq[,c("TTseq_WT_pull", "TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

colnames(merged) <- c("Gene_name","Chr","Start","End","Strand","PausIndex_WT","PausIndex_TKO",
                      "TTseq_WT_pull", "TTseq_TKO_pull", "Exon_number","Size")

# ------------------
# Save merged file |
# ------------------
write.table(merged, file = paste0(output,"TTseq_PausingIndex/TTseq_PausingIndex.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

# ------------------
# Correlation test |
# ------------------
# PausingIndex WTmean - TTseq WT pull
png(file = paste0(output, "TTseq_PausingIndex/Spearman_TTseq_PausingIndex_WT.png"))
number <- nrow(merged)
head(merged)
ggscatter(merged, y = "PausIndex_WT", x = "TTseq_WT_pull",
          color="blue",shape = 10,
          size=0.5,
          alpha=.2,
          add.params = list(color = "blue", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = F, cor.method = "spearman",
          ylim=c(0,15),
          ylab = "Pausing Index (PR/GB)", xlab = "Elongation rate (Kb/min)") +
  theme_minimal()+
  stat_cor(method = "spearman", color = "blue", label.x = 6.5, label.y = 13,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n Pausing Index vs. elongation rates in WT \n (", number, " genes)"))
dev.off()

# PausingIndex TKOmean - TTseq TKO pull
png(file = paste0(output, "TTseq_PausingIndex/Spearman_TTseq_PausingIndex_TKO.png"))
number <- nrow(merged)
head(merged)
ggscatter(merged, y = "PausIndex_TKO", x = "TTseq_TKO_pull",
          color="red2",shape = 10,
          size=0.5,
          alpha=.2,
          add.params = list(color = "red2", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = F, cor.method = "spearman",
          ylim=c(0,15),
          ylab = "Pausing Index (PR/GB)", xlab = "Elongation rate (Kb/min)") +
  theme_minimal()+
  stat_cor(method = "spearman", color = "red2", label.x = 6.5, label.y = 13,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n Pausing Index vs. elongation rates in TKO \n (", number, " genes)"))
dev.off()
