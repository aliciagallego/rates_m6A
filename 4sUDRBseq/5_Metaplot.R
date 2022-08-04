#!/usr/bin/env Rscript

#########################
## Metaplots 4sUDRB-seq #
#########################

# -------
# Paths |
# -------
WT_pull_path <- "/media/cc/B/Josemi/TTseq_LowCoverage/TTseq_LowCoverage_output/Pull_Metaplots/Histograms/WT_metaplot.txt"
TKO_pull_path <- "/media/cc/B/Josemi/TTseq_LowCoverage/TTseq_LowCoverage_output/Pull_Metaplots/Histograms/TKO_metaplot.txt"

output <- ("/media/cc/B/Josemi/TTseq_LowCoverage/TTseq_LowCoverage_output/Pull_Metaplots/Histograms/")

# -----------
# Open data |
# -----------
WT_pull <- read.table(WT_pull_path,h=T,sep="\t",stringsAsFactors=FALSE)
TKO_pull <- read.table(TKO_pull_path,h=T,sep="\t",stringsAsFactors=FALSE)

head(WT_pull)
head(TKO_pull)

tail(WT_pull)
tail(TKO_pull)

nrow(WT_pull)
nrow(TKO_pull)

# ----------
# Metaplot |
# ----------
#width=ancho
#height=alto 

png(file = paste0(output, "Metaplot_TTseq_LowCov_WT_Pull.png"))
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

png(file = paste0(output, "Metaplot_TTseq_LowCov_TKO_Pull.png"))
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


# lwd 6
plot(WT_pull$X, WT_pull$WT_15,type="l", ylim=c(0,70000),xlim=c(-2000,50000),lwd=6, col="deepskyblue3",
     xlab = "Gene body (-2Kb to +50Kb)",
     ylab = "Number of reads",
     main = "WT")
lines(WT_pull$X,WT_pull$WT_45, col="deepskyblue",lwd=6)
lines(WT_pull$X,WT_pull$WT_5, col="blue4",lwd=6)
lines(WT_pull$X,WT_pull$WT_0, col="azure4",lwd=6)
legend("topright", c("WT 0'","WT 5'","WT 15'","WT 45'"), 
       col=c("azure4","blue4","deepskyblue3","deepskyblue"), lwd=4, inset =.02, cex=0.9)

plot(TKO_pull$X, TKO_pull$TKO_15,type="l", ylim=c(0,70000),xlim=c(-2000,50000),lwd=6, col="red2",
     xlab = "Gene body (-2Kb to +50Kb)",
     ylab = "Number of reads",
     main = "H1-TKO")
lines(TKO_pull$X,TKO_pull$TKO_45, col="salmon",lwd=6)
lines(TKO_pull$X,TKO_pull$TKO_5, col="red4",lwd=6)
lines(TKO_pull$X,TKO_pull$TKO_0, col="azure4",lwd=6)

legend("topright", c("H1-TKO 0'","H1-TKO 5'","H1-TKO 15'","H1-TKO 45'"), 
       col=c("azure4","red4","red2","salmon"), lwd=4, inset =.02, cex=0.9)
