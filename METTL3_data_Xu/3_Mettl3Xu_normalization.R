#!/usr/bin/env Rscript

#########################################
## Normalize Mettl3 Xu et al. 2021 data #
#########################################
# data from:  https://doi.org/10.1038/s41586-021-03210-1

# This script processes Mettl3 data (from Xu et al. 2021) 
# and computes Mettl3 levels by quantile elongation rate catgories (slow, medium and fast)
# Boxplots are computed each elongation rate category

# -------
# Paths |
# -------
Mettl3_TSS_path <- "/path/Mettl3/Intersect/Intersect_Pull_2kbTSS.bed"
Mettl3_TTS_path <- "/path/Mettl3/Intersect/Intersect_Pull_2kbTTS.bed"
output <- "/path/Mettl3/"

# -----------
# Open data |
# -----------
Mettl3_TSS <- read.table(Mettl3_TSS_path,h=F,sep="\t",stringsAsFactors=F)
Mettl3_TTS <- read.table(Mettl3_TTS_path,h=F,sep="\t",stringsAsFactors=F)

# ----------
# Colnames |
# ----------
names_vec <- c("Chr", "Start", "End", "Gene_name", "NA1","Strand", "reads")
colnames(Mettl3_TSS) <- names_vec
colnames(Mettl3_TTS) <- names_vec

# -------
# Paths |
# -------
TTseq_path <- "/path/Rate_calculation/Elongation_rates.txt"
refseq_path <- ("/path/Input_genes_20Kb.txt")

# -----------
# Open data |
# -----------
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=FALSE)
refseq <- read.table(refseq_path,h=T,sep="\t",stringsAsFactors=FALSE)

# ----------------------------------------------------
# Transform chr20, chr21 to chrX,chrY in refseq list | 
# ----------------------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,3] == '20') {
    refseq[i,3]="X"
  }
  else if(refseq[i,3] == '21'){
    refseq[i,3]="Y"
  }
  count = count+1
}
print(count)

# -------------------------------------------------
# Remove white spaces in Gene_names in TTseq list | 
# -------------------------------------------------
TTseq$Gene_name <- gsub("\\s", "", TTseq$Gene_name)
TTseq[ TTseq == "NaN" ] <- NA

# ------------------------------
# Rates from bp/5min to Kb/min |
# ------------------------------
TTseq[,2] = TTseq[,2]/5000
TTseq[,3] = TTseq[,3]/5000

# ---------------------------------------------------
# Merge Elongation rates list and Input RefSeq list |
# ---------------------------------------------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by.x="Name", by.y="Gene_name", sort = F),
                list(refseq[,c("Chromosome","Start", "End", "Exon_number", "Orientation","Name")],
                     TTseq[,c("Gene_name","WT_Pull", "TK0_Pull")]))

# ---------------------------------------------------------
# Get transcript size in Kb based on End and Start coords |
# ---------------------------------------------------------
merged$Size <- (merged$End - merged$Start)/1000

# -------------------
# Delete NaN values |
# -------------------
merged_noNan <- na.omit(merged)  

# ------------------------------------------------
# Remove genes with elongation rate < 0.5 kb/min |
# ------------------------------------------------
merged_non05 <- merged_noNan[which(merged_noNan$WT_Pull > 0.5 & merged_noNan$TK0_Pull > 0.5),]
colnames(merged_non05)[1] <- "Gene_name"

# ----------------------------
# Merge rates and other data |
# ----------------------------
merged_TSS<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                           list(Mettl3_TSS[,c("Chr","Start","End","Strand","reads","Gene_name")],
                                merged_non05[,c("WT_Pull", "Exon_number", "Size","Gene_name")]))

merged_TTS<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                    list(Mettl3_TTS[,c("Chr","Start","End","Strand","reads","Gene_name")],
                         merged_non05[,c("WT_Pull", "Exon_number", "Size","Gene_name")]))

# ------------------------------
# Obtain tertiles and classify |
# ------------------------------

# WT
QWT <- quantile(merged_non05$WT_Pull, probs = seq(0, 1, 1/3), na.rm = FALSE)  
upWT <-  QWT[[3]]  
lowWT<- QWT[[2]]

WTslow <- (subset(merged_non05, merged_non05$WT_Pull<lowWT))
WTslow <- WTslow[order(WTslow$WT_Pull),]
WTmedium  <- (subset(merged_non05, merged_non05$WT_Pull>=lowWT & merged_non05$WT_Pull<=upWT))
WTmedium <- WTmedium[order(WTmedium$WT_Pull),]
WTfast <- (subset(merged_non05, merged_non05$WT_Pull>upWT))
WTfast <- WTfast[order(WTfast$WT_Pull),]

# TKO
QTKO <- quantile(merged_non05$TK0_Pull, probs = seq(0, 1, 1/3), na.rm = FALSE)  
upTKO <-  QTKO[[3]]  
lowTKO<- QTKO[[2]]

TKOslow <- (subset(merged_non05, merged_non05$TK0_Pull<lowTKO))
TKOslow <- TKOslow[order(TKOslow$TK0_Pull),]
TKOmedium  <- (subset(merged_non05, merged_non05$TK0_Pull>=lowTKO & merged_non05$TK0_Pull<=upTKO))
TKOmedium <- TKOmedium[order(TKOmedium$TK0_Pull),]
TKOfast <- (subset(merged_non05, merged_non05$TK0_Pull>upTKO))
TKOfast <- TKOfast[order(TKOfast$TK0_Pull),]

# --------------------------------
# Merged tertiles and other data |
# --------------------------------
names_vec <- c("Gene_name", "Chr", "Start", "End", "Exon_number","Strand", "TTseq_WT_pull", 
               "TTseq_TKO_pull","Size")

colnames(WTslow) <- names_vec
colnames(WTmedium) <- names_vec
colnames(WTfast) <- names_vec
colnames(TKOslow) <- names_vec
colnames(TKOmedium) <- names_vec
colnames(TKOfast) <- names_vec

merged_WTslow_TSS<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                 list(Mettl3_TSS[,c("Chr","Start","End","Strand","reads","Gene_name")],
                                      WTslow[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))
merged_WTmed_TSS<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(Mettl3_TSS[,c("Chr","Start","End","Strand","reads","Gene_name")],
                                WTmedium[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))
merged_WTfast_TSS<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                          list(Mettl3_TSS[,c("Chr","Start","End","Strand","reads","Gene_name")],
                               WTfast[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))
merged_WTslow_TTS<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(Mettl3_TTS[,c("Chr","Start","End","Strand","reads","Gene_name")],
                                WTslow[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))
merged_WTmed_TTS<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                          list(Mettl3_TTS[,c("Chr","Start","End","Strand","reads","Gene_name")],
                               WTmedium[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))
merged_WTfast_TTS<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(Mettl3_TTS[,c("Chr","Start","End","Strand","reads","Gene_name")],
                                WTfast[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

# ----------
# Boxplots |
# ----------
boxplot(merged_WTslow_TSS$reads, merged_WTmed_TSS$reads,merged_WTfast_TSS$reads,
        merged_WTslow_TTS$reads, merged_WTmed_TTS$reads,merged_WTfast_TTS$reads,
        col = "grey",
        at = c(1,2,3,5,6,7),
        names = c("TSS slow","TSS med", "TSS fast", 
                  "TTS slow","TTS med", "TTS fast"),
        main = "Mettl3 levels ( data from Xu et al.)",
        ylab="Number of reads",
        outline = F,
        las=2,
        cex=.8,
        ylim=c(-10,350),
        boxwex = 0.7,
        lwd=1.5)
segments(1,260,2,260,lwd=2)
segments(1,300,3,300,lwd=2)
text(1.5, 270,"p-val = 1.676e-04",cex=.8, col=("black"))
text(2, 310, "p-val = 2.032e-08",cex=.8,col=("black"))
dev.off()	
