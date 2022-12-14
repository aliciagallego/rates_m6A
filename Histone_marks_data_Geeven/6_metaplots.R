#!/usr/bin/env Rscript

##########################################
## Metaplots histones Geeven et al. 2015 #
##########################################

# This script generates meta-gene profiles from histograms
# Meta-gene profiles include information from H3K4me3, H3K9me3, H3K27me3 histone marks (Geeven et al. 2015) and 
# from MNase-seq (own data)

# -------
# Paths |
# -------
H3K4me3_path <- "/path/Histograms/H3K4me3/"
H3K9me3_path <- "/path/Histograms/H3K9me3/"
H3K27me3_path <- "/path/Histograms/H3K27me3/"

WT_MNase_path <- "/media/cc/C/Alicia/MNase/histograms/WT_metaplot.txt"
TKO_MNase_path <- "/media/cc/C/Alicia/MNase/histograms/TKO_metaplot.txt"

output <- ("/path/Histograms/")

# --------------------------------------------------------
# Read all files from directories and put them in a list |
# --------------------------------------------------------
WT_MNase <- read.table(WT_MNase_path,h=T,sep="\t",stringsAsFactors=F)
TKO_MNase <- read.table(TKO_MNase_path,h=T,sep="\t",stringsAsFactors=F)

H3K4me3_list = list.files(H3K4me3_path, pattern="*txt")
for (i in seq_along(H3K4me3_list)) {
  filename <- sub(".txt", "", H3K4me3_list[i])
  df <- read.table(paste0(H3K4me3_path,H3K4me3_list[i]), header = T, sep="\t", stringsAsFactors=FALSE)
  assign(filename, df)
}

H3K9me3_list = list.files(H3K9me3_path, pattern="*txt")
for (i in seq_along(H3K9me3_list)) {
  filename <- sub(".txt", "", H3K9me3_list[i])
  df <- read.table(paste0(H3K9me3_path,H3K9me3_list[i]), header = T, sep="\t", stringsAsFactors=FALSE)
  assign(filename, df)
}

H3K27me3_list = list.files(H3K27me3_path, pattern="*txt")
for (i in seq_along(H3K27me3_list)) {
  filename <- sub(".txt", "", H3K27me3_list[i])
  df <- read.table(paste0(H3K27me3_path,H3K27me3_list[i]), header = T, sep="\t", stringsAsFactors=FALSE)
  assign(filename, df)
}

# ----------------
# Transform data | relative to each mark
# ----------------
H3K4me3_WT_total <- sum(H3K4me3_WT$slow)+sum(H3K4me3_WT$med)+sum(H3K4me3_WT$fast)
H3K4me3_WT$slow <- H3K4me3_WT$slow/H3K4me3_WT_total*100
H3K4me3_WT$med <- H3K4me3_WT$med/H3K4me3_WT_total*100
H3K4me3_WT$fast <- H3K4me3_WT$fast/H3K4me3_WT_total*100

H3K4me3_TKO_total <- sum(H3K4me3_TKO$slow)+sum(H3K4me3_TKO$med)+sum(H3K4me3_TKO$fast)
H3K4me3_TKO$slow <- H3K4me3_TKO$slow/H3K4me3_TKO_total*100
H3K4me3_TKO$med <- H3K4me3_TKO$med/H3K4me3_TKO_total*100
H3K4me3_TKO$fast <- H3K4me3_TKO$fast/H3K4me3_TKO_total*100

H3K9me3_WT_total <- sum(H3K9me3_WT$slow)+sum(H3K9me3_WT$med)+sum(H3K9me3_WT$fast)
H3K9me3_WT$slow <- H3K9me3_WT$slow/H3K9me3_WT_total*100
H3K9me3_WT$med <- H3K9me3_WT$med/H3K9me3_WT_total*100
H3K9me3_WT$fast <- H3K9me3_WT$fast/H3K9me3_WT_total*100

H3K9me3_TKO_total <- sum(H3K9me3_TKO$slow)+sum(H3K9me3_TKO$med)+sum(H3K9me3_TKO$fast)
H3K9me3_TKO$slow <- H3K9me3_TKO$slow/H3K9me3_TKO_total*100
H3K9me3_TKO$med <- H3K9me3_TKO$med/H3K9me3_TKO_total*100
H3K9me3_TKO$fast <- H3K9me3_TKO$fast/H3K9me3_TKO_total*100

H3K27me3_WT_total <- sum(H3K27me3_WT$slow)+sum(H3K27me3_WT$med)+sum(H3K27me3_WT$fast)
H3K27me3_WT$slow <- H3K27me3_WT$slow/H3K27me3_WT_total*100
H3K27me3_WT$med <- H3K27me3_WT$med/H3K27me3_WT_total*100
H3K27me3_WT$fast <- H3K27me3_WT$fast/H3K27me3_WT_total*100

H3K27me3_TKO_total <- sum(H3K27me3_TKO$slow)+sum(H3K27me3_TKO$med)+sum(H3K27me3_TKO$fast)
H3K27me3_TKO$slow <- H3K27me3_TKO$slow/H3K27me3_TKO_total*100
H3K27me3_TKO$med <- H3K27me3_TKO$med/H3K27me3_TKO_total*100
H3K27me3_TKO$fast <- H3K27me3_TKO$fast/H3K27me3_TKO_total*100

WT_MNase_total <- sum(WT_MNase$WT_slow)+sum(WT_MNase$WT_medium)+sum(WT_MNase$WT_fast)
WT_MNase$WT_slow <- WT_MNase$WT_slow/WT_MNase_total*100
WT_MNase$WT_medium <- WT_MNase$WT_medium/WT_MNase_total*100
WT_MNase$WT_fast <- WT_MNase$WT_fast/WT_MNase_total*100

TKO_MNase_total <- sum(TKO_MNase$TKO_slow)+sum(TKO_MNase$TKO_medium)+sum(TKO_MNase$TKO_fast)
TKO_MNase$TKO_slow <- TKO_MNase$TKO_slow/TKO_MNase_total*100
TKO_MNase$TKO_medium <- TKO_MNase$TKO_medium/TKO_MNase_total*100
TKO_MNase$TKO_fast <- TKO_MNase$TKO_fast/TKO_MNase_total*100

# --------------------
# Meta-gene profiles |
# --------------------

# WT slow
plot(H3K4me3_WT$X, H3K4me3_WT$slow,type="l", ylim=c(0,0.6), lwd=6, col="darkgreen",
     xlab = "TSS??2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications and nucleosome profiles in WT slow rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K9me3_WT$X,H3K9me3_WT$slow, col="dodgerblue3",lwd=6)
lines(H3K27me3_WT$X,H3K27me3_WT$slow, col="darkslateblue",lwd=6)
lines(WT_MNase$X, WT_MNase$WT_slow, col="darkolivegreen3",lwd=6)
legend("topright", c("H3K4me3 WT","H3K9me3 WT","H3K27me3 WT","MNase WT"), 
       col=c("darkgreen","dodgerblue3","darkslateblue","darkolivegreen3"), 
       lwd=4, inset =.01, cex=1, box.col = "white")

# WT med
plot(H3K4me3_WT$X, H3K4me3_WT$med,type="l", ylim=c(0,0.6), lwd=6, col="darkgreen",
     xlab = "TSS??2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications and nucleosome profiles in WT medium rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K9me3_WT$X,H3K9me3_WT$med, col="dodgerblue3",lwd=6)
lines(H3K27me3_WT$X,H3K27me3_WT$med, col="darkslateblue",lwd=6)
lines(WT_MNase$X, WT_MNase$WT_medium, col="darkolivegreen3",lwd=6)
legend("topright", c("H3K4me3 WT","H3K9me3 WT","H3K27me3 WT","MNase WT"), 
       col=c("darkgreen","dodgerblue3","darkslateblue","darkolivegreen3"), 
       lwd=4, inset =.01, cex=1,box.col = "white")

# WT fast
plot(H3K4me3_WT$X, H3K4me3_WT$fast,type="l", ylim=c(0,0.6), lwd=6, col="darkgreen",
     xlab = "TSS??2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications and nucleosome profiles in WT fast rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K9me3_WT$X,H3K9me3_WT$fast, col="dodgerblue3",lwd=6)
lines(H3K27me3_WT$X,H3K27me3_WT$fast, col="darkslateblue",lwd=6)
lines(WT_MNase$X, WT_MNase$WT_fast, col="darkolivegreen3",lwd=6)
legend("topright", c("H3K4me3 WT","H3K9me3 WT","H3K27me3 WT","MNase WT"), 
       col=c("darkgreen","dodgerblue3","darkslateblue","darkolivegreen3"), 
       lwd=4, inset =.01, cex=1,box.col = "white")

# TKO slow
plot(H3K4me3_TKO$X, H3K4me3_TKO$slow,type="l", ylim=c(0,0.6), lwd=6, col="darkorange2",
     xlab = "TSS??2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications and nucleosome profiles in TKO slow rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K9me3_TKO$X,H3K9me3_TKO$slow, col="firebrick3",lwd=6)
lines(H3K27me3_TKO$X,H3K27me3_TKO$slow, col="deeppink3",lwd=6)
lines(TKO_MNase$X, TKO_MNase$TKO_slow, col="darkgoldenrod1",lwd=6)
legend("topright", c("H3K4me3 TKO","H3K9me3 TKO","H3K27me3 TKO","MNase TKO"), 
       col=c("darkorange2","firebrick3","deeppink3","darkgoldenrod1"), 
       lwd=4, inset =.01, cex=1,box.col = "white")

# TKO med
plot(H3K4me3_TKO$X, H3K4me3_TKO$med,type="l", ylim=c(0,0.6), lwd=6, col="darkorange2",
     xlab = "TSS??2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications and nucleosome profiles in TKO medium rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K9me3_TKO$X,H3K9me3_TKO$med, col="firebrick3",lwd=6)
lines(H3K27me3_TKO$X,H3K27me3_TKO$med, col="deeppink3",lwd=6)
lines(TKO_MNase$X, TKO_MNase$TKO_medium, col="darkgoldenrod1",lwd=6)
legend("topright", c("H3K4me3 TKO","H3K9me3 TKO","H3K27me3 TKO","MNase TKO"), 
       col=c("darkorange2","firebrick3","deeppink3","darkgoldenrod1"), 
       lwd=4, inset =.01, cex=1,box.col = "white")

# TKO fast
plot(H3K4me3_TKO$X, H3K4me3_TKO$fast,type="l", ylim=c(0,0.6), lwd=6, col="darkorange2",
     xlab = "TSS??2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications and nucleosome profiles in TKO fast rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K9me3_TKO$X,H3K9me3_TKO$fast, col="firebrick3",lwd=6)
lines(H3K27me3_TKO$X,H3K27me3_TKO$fast, col="deeppink3",lwd=6)
lines(TKO_MNase$X, TKO_MNase$TKO_fast, col="darkgoldenrod1",lwd=6)
legend("topright", c("H3K4me3 TKO","H3K9me3 TKO","H3K27me3 TKO","MNase TKO"), 
       col=c("darkorange2","firebrick3","deeppink3","darkgoldenrod1"), 
       lwd=4, inset =.01, cex=1,box.col = "white")
