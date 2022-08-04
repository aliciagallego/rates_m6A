#!/usr/bin/env Rscript

##########################
## Metaplots 4sU-DRB-seq #
##########################

# This script generates meta-gene profiles from histograms
# Histogram information from 0, 5, 15 and 45 min (for both WT and H1-TKO samples) were previously pasted together in the same file

# -------
# Paths |
# -------
WT_path <- "/path/Histograms/WT_Pull_matrix_sum_001.txt"
TKO_path <- "/path/Histograms/TKO_Pull_matrix_sum_001.txt.txt"
output <- "path/Histograms/"

# -----------
# Open data |
# -----------
WT_pull <- read.table(WT_pull_path,h=T,sep="\t",stringsAsFactors=FALSE)
TKO_pull <- read.table(TKO_pull_path,h=T,sep="\t",stringsAsFactors=FALSE)

# -----------------
# Meta-gene plots |
# -----------------

# WT profiles
png(file = paste0(output, "Metaplot_TTseq_WT_Pull.png"))
plot(WT_pull$X, WT_pull$WT_5,type="l", ylim=c(0,70000),xlim=c(-2000,50000),lwd=3, col="blue4",
     xlab = "Gene body (-2Kb to +50Kb)",
     ylab = "Number of reads",
     main = "WT")
lines(WT_pull$X,WT_pull$WT_0, col="grey",lwd=3)
lines(WT_pull$X,WT_pull$WT_15, col="deepskyblue3",lwd=3)
lines(WT_pull$X,WT_pull$WT_45, col="deepskyblue",lwd=3)
legend("topright", c("WT 0'","WT 5'","WT 15'","WT 45'"), 
       col=c("grey","blue4","deepskyblue3","deepskyblue"), lwd=3, inset =.02, cex=0.9)
dev.off()

# H1-TKO profiles
png(file = paste0(output, "Metaplot_TTseq_TKO_Pull.png"))
plot(TKO_pull$X, TKO_pull$TKO_5,type="l", ylim=c(0,70000),xlim=c(-2000,50000),lwd=3, col="red4",
     xlab = "Gene body (-2Kb to +50Kb)",
     ylab = "Number of reads",
     main = "H1-TKO")
lines(TKO_pull$X,TKO_pull$TKO_0, col="grey",lwd=3)
lines(TKO_pull$X,TKO_pull$TKO_15, col="red2",lwd=3)
lines(TKO_pull$X,TKO_pull$TKO_45, col="salmon",lwd=3)
legend("topright", c("H1-TKO 0'","H1-TKO 5'","H1-TKO 15'","H1-TKO 45'"), 
       col=c("grey","red4","red2","salmon"), lwd=3, inset =.02, cex=0.9)
dev.off()
