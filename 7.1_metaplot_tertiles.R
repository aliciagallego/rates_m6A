#!/usr/bin/env Rscript

##########################################
## Metaplots histones Geeven et al. 2015 #
##########################################

# -------
# Paths |
# -------
H3K4me1_path <- "/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/H3K4me1/"
H3K4me3_path <- "/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/H3K4me3/"
H3K9me3_path <- "/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/H3K9me3/"
H3K27me3_path <- "/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/H3K27me3/"

WT_MNase_path <- "/media/cc/C/Alicia/MNase/histograms/WT_metaplot.txt"
TKO_MNase_path <- "/media/cc/C/Alicia/MNase/histograms/TKO_metaplot.txt"

output <- ("/media/cc/B/Alicia/Geeven/Geeven_output/3.2_Histograms_tertiles/")

# --------------------------------------------------------
# Read all files from directories and put them in a list |
# --------------------------------------------------------
WT_MNase <- read.table(WT_MNase_path,h=T,sep="\t",stringsAsFactors=F)
TKO_MNase <- read.table(TKO_MNase_path,h=T,sep="\t",stringsAsFactors=F)

H3K4me1_list = list.files(H3K4me1_path, pattern="*txt")
for (i in seq_along(H3K4me1_list)) {
  filename <- sub(".txt", "", H3K4me1_list[i])
  df <- read.table(paste0(H3K4me1_path,H3K4me1_list[i]), header = T, sep="\t", stringsAsFactors=FALSE)
  assign(filename, df)
}

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
head(H3K4me1_WT)
H3K4me1_WT_total <- sum(H3K4me1_WT$slow)+sum(H3K4me1_WT$med)+sum(H3K4me1_WT$fast)
H3K4me1_WT_total
H3K4me1_WT$slow <- H3K4me1_WT$slow/H3K4me1_WT_total*100
H3K4me1_WT$med <- H3K4me1_WT$med/H3K4me1_WT_total*100
H3K4me1_WT$fast <- H3K4me1_WT$fast/H3K4me1_WT_total*100

H3K4me1_TKO_total <- sum(H3K4me1_TKO$slow)+sum(H3K4me1_TKO$med)+sum(H3K4me1_TKO$fast)
H3K4me1_TKO_total
H3K4me1_TKO$slow <- H3K4me1_TKO$slow/H3K4me1_TKO_total*100
H3K4me1_TKO$med <- H3K4me1_TKO$med/H3K4me1_TKO_total*100
H3K4me1_TKO$fast <- H3K4me1_TKO$fast/H3K4me1_TKO_total*100

H3K4me3_WT_total <- sum(H3K4me3_WT$slow)+sum(H3K4me3_WT$med)+sum(H3K4me3_WT$fast)
H3K4me3_WT_total
H3K4me3_WT$slow <- H3K4me3_WT$slow/H3K4me3_WT_total*100
H3K4me3_WT$med <- H3K4me3_WT$med/H3K4me3_WT_total*100
H3K4me3_WT$fast <- H3K4me3_WT$fast/H3K4me3_WT_total*100

H3K4me3_TKO_total <- sum(H3K4me3_TKO$slow)+sum(H3K4me3_TKO$med)+sum(H3K4me3_TKO$fast)
H3K4me3_TKO_total
H3K4me3_TKO$slow <- H3K4me3_TKO$slow/H3K4me3_TKO_total*100
H3K4me3_TKO$med <- H3K4me3_TKO$med/H3K4me3_TKO_total*100
H3K4me3_TKO$fast <- H3K4me3_TKO$fast/H3K4me3_TKO_total*100

H3K9me3_WT_total <- sum(H3K9me3_WT$slow)+sum(H3K9me3_WT$med)+sum(H3K9me3_WT$fast)
H3K9me3_WT_total
H3K9me3_WT$slow <- H3K9me3_WT$slow/H3K9me3_WT_total*100
H3K9me3_WT$med <- H3K9me3_WT$med/H3K9me3_WT_total*100
H3K9me3_WT$fast <- H3K9me3_WT$fast/H3K9me3_WT_total*100

H3K9me3_TKO_total <- sum(H3K9me3_TKO$slow)+sum(H3K9me3_TKO$med)+sum(H3K9me3_TKO$fast)
H3K9me3_TKO_total
H3K9me3_TKO$slow <- H3K9me3_TKO$slow/H3K9me3_TKO_total*100
H3K9me3_TKO$med <- H3K9me3_TKO$med/H3K9me3_TKO_total*100
H3K9me3_TKO$fast <- H3K9me3_TKO$fast/H3K9me3_TKO_total*100

H3K27me3_WT_total <- sum(H3K27me3_WT$slow)+sum(H3K27me3_WT$med)+sum(H3K27me3_WT$fast)
H3K27me3_WT_total
H3K27me3_WT$slow <- H3K27me3_WT$slow/H3K27me3_WT_total*100
H3K27me3_WT$med <- H3K27me3_WT$med/H3K27me3_WT_total*100
H3K27me3_WT$fast <- H3K27me3_WT$fast/H3K27me3_WT_total*100

H3K27me3_TKO_total <- sum(H3K27me3_TKO$slow)+sum(H3K27me3_TKO$med)+sum(H3K27me3_TKO$fast)
H3K27me3_TKO_total
H3K27me3_TKO$slow <- H3K27me3_TKO$slow/H3K27me3_TKO_total*100
H3K27me3_TKO$med <- H3K27me3_TKO$med/H3K27me3_TKO_total*100
H3K27me3_TKO$fast <- H3K27me3_TKO$fast/H3K27me3_TKO_total*100

head(WT_MNase)
WT_MNase_total <- sum(WT_MNase$WT_slow)+sum(WT_MNase$WT_medium)+sum(WT_MNase$WT_fast)
WT_MNase_total
WT_MNase$WT_slow <- WT_MNase$WT_slow/WT_MNase_total*100
WT_MNase$WT_medium <- WT_MNase$WT_medium/WT_MNase_total*100
WT_MNase$WT_fast <- WT_MNase$WT_fast/WT_MNase_total*100

head(TKO_MNase)
TKO_MNase_total <- sum(TKO_MNase$TKO_slow)+sum(TKO_MNase$TKO_medium)+sum(TKO_MNase$TKO_fast)
TKO_MNase_total
TKO_MNase$TKO_slow <- TKO_MNase$TKO_slow/TKO_MNase_total*100
TKO_MNase$TKO_medium <- TKO_MNase$TKO_medium/TKO_MNase_total*100
TKO_MNase$TKO_fast <- TKO_MNase$TKO_fast/TKO_MNase_total*100

# ----------
# Metaplot |
# ----------
#width=ancho
#height=alto 

par(mar=c(4.2, 4.2, 4, 2))
# con H3K4me1
# WT slow
plot(H3K4me1_WT$X, H3K4me1_WT$slow,type="l", ylim=c(0,0.5), lwd=6, col="darkolivegreen3",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications around TSS in WT slow rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K4me3_WT$X,H3K4me3_WT$slow, col="darkgreen",lwd=4)
lines(H3K9me3_WT$X,H3K9me3_WT$slow, col="dodgerblue3",lwd=4)
lines(H3K27me3_WT$X,H3K27me3_WT$slow, col="darkslateblue",lwd=4)
legend("topright", c("H3K4me1 WT","H3K4me3 WT","H3K9me3 WT","H3K27me3 WT"), 
       col=c("darkolivegreen3","darkgreen","dodgerblue3","darkslateblue"), 
       lwd=3, inset =.01, cex=1, box.col = "white")

# WT med
plot(H3K4me1_WT$X, H3K4me1_WT$med,type="l", ylim=c(0,0.5), lwd=4, col="darkolivegreen3",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications around TSS in WT med rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K4me3_WT$X,H3K4me3_WT$med, col="darkgreen",lwd=4)
lines(H3K9me3_WT$X,H3K9me3_WT$med, col="dodgerblue3",lwd=4)
lines(H3K27me3_WT$X,H3K27me3_WT$med, col="darkslateblue",lwd=4)
legend("topright", c("H3K4me1 WT","H3K4me3 WT","H3K9me3 WT","H3K27me3 WT"), 
       col=c("darkolivegreen3","darkgreen","dodgerblue3","darkslateblue"), 
       lwd=3, inset =.01, cex=1,box.col = "white")

# WT fast
plot(H3K4me1_WT$X, H3K4me1_WT$fast,type="l", ylim=c(0,0.5), lwd=4, col="darkolivegreen3",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications around TSS in WT fast rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K4me3_WT$X,H3K4me3_WT$fast, col="darkgreen",lwd=4)
lines(H3K9me3_WT$X,H3K9me3_WT$fast, col="dodgerblue3",lwd=4)
lines(H3K27me3_WT$X,H3K27me3_WT$fast, col="darkslateblue",lwd=4)
legend("topright", c("H3K4me1 WT","H3K4me3 WT","H3K9me3 WT","H3K27me3 WT"), 
       col=c("darkolivegreen3","darkgreen","dodgerblue3","darkslateblue"), 
       lwd=3, inset =.01, cex=1,box.col = "white")


# TKO slow
plot(H3K4me1_TKO$X, H3K4me1_TKO$slow,type="l", ylim=c(0,0.5), lwd=4, col="darkgoldenrod1",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications around TSS in TKO slow rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K4me3_TKO$X,H3K4me3_TKO$slow, col="darkorange2",lwd=4)
lines(H3K9me3_TKO$X,H3K9me3_TKO$slow, col="firebrick3",lwd=4)
lines(H3K27me3_TKO$X,H3K27me3_TKO$slow, col="deeppink3",lwd=4)
legend("topright", c("H3K4me1 TKO","H3K4me3 TKO","H3K9me3 TKO","H3K27me3 TKO"), 
       col=c("darkgoldenrod1","darkorange2","firebrick3","deeppink3"), 
       lwd=3, inset =.01, cex=1,box.col = "white")

# TKO med
plot(H3K4me1_TKO$X, H3K4me1_TKO$med,type="l", ylim=c(0,0.5), lwd=4, col="darkgoldenrod1",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications around TSS in TKO med rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K4me3_TKO$X,H3K4me3_TKO$med, col="darkorange2",lwd=4)
lines(H3K9me3_TKO$X,H3K9me3_TKO$med, col="firebrick3",lwd=4)
lines(H3K27me3_TKO$X,H3K27me3_TKO$med, col="deeppink3",lwd=4)
legend("topright", c("H3K4me1 TKO","H3K4me3 TKO","H3K9me3 TKO","H3K27me3 TKO"), 
       col=c("darkgoldenrod1","darkorange2","firebrick3","deeppink3"), 
       lwd=3, inset =.01, cex=1,box.col = "white")

# TKO fast
plot(H3K4me1_TKO$X, H3K4me1_TKO$fast,type="l", ylim=c(0,0.5), lwd=4, col="darkgoldenrod1",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications around TSS in TKO fast rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K4me3_TKO$X,H3K4me3_TKO$fast, col="darkorange2",lwd=4)
lines(H3K9me3_TKO$X,H3K9me3_TKO$fast, col="firebrick3",lwd=4)
lines(H3K27me3_TKO$X,H3K27me3_TKO$fast, col="deeppink3",lwd=4)
legend("topright", c("H3K4me1 TKO","H3K4me3 TKO","H3K9me3 TKO","H3K27me3 TKO"), 
       col=c("darkgoldenrod1","darkorange2","firebrick3","deeppink3"), 
       lwd=3, inset =.01, cex=1,box.col = "white")

# sin H3K4me1
# WT slow
plot(H3K4me3_WT$X, H3K4me3_WT$slow,type="l", ylim=c(0,0.6), lwd=6, col="darkgreen",
     xlab = "TSS±2Kb",
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
     xlab = "TSS±2Kb",
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
     xlab = "TSS±2Kb",
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
     xlab = "TSS±2Kb",
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
     xlab = "TSS±2Kb",
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
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Histone modifications and nucleosome profiles in TKO fast rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K9me3_TKO$X,H3K9me3_TKO$fast, col="firebrick3",lwd=6)
lines(H3K27me3_TKO$X,H3K27me3_TKO$fast, col="deeppink3",lwd=6)
lines(TKO_MNase$X, TKO_MNase$TKO_fast, col="darkgoldenrod1",lwd=6)
legend("topright", c("H3K4me3 TKO","H3K9me3 TKO","H3K27me3 TKO","MNase TKO"), 
       col=c("darkorange2","firebrick3","deeppink3","darkgoldenrod1"), 
       lwd=4, inset =.01, cex=1,box.col = "white")

# MNase

plot(WT_MNase$X, WT_MNase$WT_fast,type="l", ylim=c(0,0.6), lwd=6, col="blue4",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Nucleosome profiles in WT",
     cex.lab = 1.5,cex.axis = 1.5)
lines(WT_MNase$X,WT_MNase$WT_medium, col="deepskyblue3",lwd=6)
lines(WT_MNase$X,WT_MNase$WT_slow, col="lightblue3",lwd=6)
legend("topright", c("WT slow","WT medium","WT fast"), 
       col=c("lightblue3","deepskyblue3","blue4"), 
       lwd=4, inset =.01, cex=1,box.col = "white")

plot(TKO_MNase$X, TKO_MNase$TKO_fast,type="l", ylim=c(0,0.6), lwd=6, col="red3",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Nucleosome profiles in WT",
     cex.lab = 1.5,cex.axis = 1.5)
lines(TKO_MNase$X,TKO_MNase$TKO_medium, col="red",lwd=6)
lines(TKO_MNase$X,TKO_MNase$TKO_slow, col="salmon",lwd=6)
legend("topright", c("TKO slow","TKO medium","TKO fast"), 
       col=c("salmon","red","red3"),
       lwd=4, inset =.01, cex=1,box.col = "white")

plot(WT_MNase$X, WT_MNase$WT_slow,type="l", ylim=c(0,0.6), lwd=6, col="lightblue3",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Nucleosome profiles in slow rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(TKO_MNase$X,TKO_MNase$TKO_slow, col="salmon",lwd=6)
legend("topright", c("WT slow","TKO slow"), 
       col=c("lightblue3","salmon"), 
       lwd=4, inset =.01, cex=1,box.col = "white")

plot(WT_MNase$X, WT_MNase$WT_medium,type="l", ylim=c(0,0.6), lwd=6, col="deepskyblue3",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Nucleosome profiles in medium rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(TKO_MNase$X,TKO_MNase$TKO_medium, col="red",lwd=6)
legend("topright", c("WT slow","TKO slow"), 
       col=c("deepskyblue3","red"), 
       lwd=4, inset =.01, cex=1,box.col = "white")

plot(WT_MNase$X, WT_MNase$WT_fast,type="l", ylim=c(0,0.6), lwd=6, col="blue4",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "Nucleosome profiles in fast rate genes",
     cex.lab = 1.5,cex.axis = 1.5)
lines(TKO_MNase$X,TKO_MNase$TKO_fast, col="red3",lwd=6)
legend("topright", c("WT slow","TKO slow"), 
       col=c("blue4","red3"), 
       lwd=4, inset =.01, cex=1,box.col = "white")

# By histone marks WT
# H3K4me3
plot(H3K4me3_WT$X, H3K4me3_WT$slow,type="l", ylim=c(0,0.6), lwd=6, col="deepskyblue",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "H3K4me3 histone mark in WT",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K4me3_WT$X,H3K4me3_WT$med, col="deepskyblue3",lwd=6)
lines(H3K4me3_WT$X,H3K4me3_WT$fast, col="dodgerblue4",lwd=6)
legend("topright", c("WT slow","WT medium","WT fast"), 
       col=c("deepskyblue","deepskyblue3","dodgerblue4"), 
       lwd=4, inset =.01, cex=1, box.col = "white")

# H3K9me3
plot(H3K9me3_WT$X, H3K9me3_WT$slow,type="l", ylim=c(0,0.6), lwd=6, col="deepskyblue",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "H3K9me3 histone mark in WT",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K9me3_WT$X,H3K9me3_WT$med, col="deepskyblue3",lwd=6)
lines(H3K9me3_WT$X,H3K9me3_WT$fast, col="dodgerblue4",lwd=6)
legend("topright", c("WT slow","WT medium","WT fast"), 
       col=c("deepskyblue","deepskyblue3","dodgerblue4"), 
       lwd=4, inset =.01, cex=1, box.col = "white")

# H3K27me3
plot(H3K27me3_WT$X, H3K27me3_WT$slow,type="l", ylim=c(0,0.6), lwd=6, col="deepskyblue",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "H3K27me3 histone mark in WT",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K27me3_WT$X,H3K27me3_WT$med, col="deepskyblue3",lwd=6)
lines(H3K27me3_WT$X,H3K27me3_WT$fast, col="dodgerblue4",lwd=6)
legend("topright", c("WT slow","WT medium","WT fast"), 
       col=c("deepskyblue","deepskyblue3","dodgerblue4"), 
       lwd=4, inset =.01, cex=1, box.col = "white")

# MNase
plot(WT_MNase$X, WT_MNase$WT_slow,type="l", ylim=c(0,0.6), lwd=6, col="deepskyblue",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "MNase histone mark in WT",
     cex.lab = 1.5,cex.axis = 1.5)
lines(WT_MNase$X,WT_MNase$WT_medium, col="deepskyblue3",lwd=6)
lines(WT_MNase$X,WT_MNase$WT_fast, col="dodgerblue4",lwd=6)
legend("topright", c("WT slow","WT medium","WT fast"), 
       col=c("deepskyblue","deepskyblue3","dodgerblue4"), 
       lwd=4, inset =.01, cex=1, box.col = "white")

# By histone marks TKO
# H3K4me3
plot(H3K4me3_TKO$X, H3K4me3_TKO$slow,type="l", ylim=c(0,0.6), lwd=6, col="salmon",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "H3K4me3 histone mark in TKO",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K4me3_TKO$X,H3K4me3_TKO$med, col="darkorange",lwd=6)
lines(H3K4me3_TKO$X,H3K4me3_TKO$fast, col="red",lwd=6)
legend("topright", c("TKO slow","TKO medium","TKO fast"), 
       col=c("salmon","darkorange","red"), 
       lwd=4, inset =.01, cex=1, box.col = "white")

# H3K9me3
plot(H3K9me3_TKO$X, H3K9me3_TKO$slow,type="l", ylim=c(0,0.6), lwd=6, col="salmon",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "H3K9me3 histone mark in TKO",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K9me3_TKO$X,H3K9me3_TKO$med, col="darkorange",lwd=6)
lines(H3K9me3_TKO$X,H3K9me3_TKO$fast, col="red",lwd=6)
legend("topright", c("TKO slow","TKO medium","TKO fast"), 
       col=c("salmon","darkorange","red"), 
       lwd=4, inset =.01, cex=1, box.col = "white")

# H3K27me3
plot(H3K27me3_TKO$X, H3K27me3_TKO$slow,type="l", ylim=c(0,0.6), lwd=6, col="salmon",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "H3K27me3 histone mark in TKO",
     cex.lab = 1.5,cex.axis = 1.5)
lines(H3K27me3_TKO$X,H3K27me3_TKO$med, col="darkorange",lwd=6)
lines(H3K27me3_TKO$X,H3K27me3_TKO$fast, col="red",lwd=6)
legend("topright", c("TKO slow","TKO medium","TKO fast"), 
       col=c("salmon","darkorange","red"), 
       lwd=4, inset =.01, cex=1, box.col = "white")

# MNase
plot(TKO_MNase$X, TKO_MNase$TKO_slow,type="l", xlim=c(-2000,2000), ylim=c(0,0.6), lwd=6, col="salmon",
     xlab = "TSS±2Kb",
     ylab = "Normalized reads",
     main = "MNase histone mark in TKO",
     cex.lab = 1.5,cex.axis = 1.5)
lines(TKO_MNase$X,TKO_MNase$TKO_medium, col="darkorange",lwd=6)
lines(TKO_MNase$X,TKO_MNase$TKO_fast, col="red",lwd=6)
legend("topright", c("TKO slow","TKO medium","TKO fast"), 
       col=c("salmon","darkorange","red"), 
       lwd=4, inset =.01, cex=1, box.col = "white")
