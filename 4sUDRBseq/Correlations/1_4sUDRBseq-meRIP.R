#!/usr/bin/env Rscript

##########################################################################
## Correlations elongation rates (4sUDRBseq or TTseq) and m6A (meRIPseq) #
##########################################################################

# -----------
# Libraries |
# -----------
library("ggpubr")
library("ggplot2")
library("reshape2")

# -------
# Paths |
# -------
m6AWT_path <- "/path/meRIP/Normalized_data/meRIP_RefSeqLongList_Normalized_WT_2Kb.txt"
m6ATKO_path <-"/path/meRIP/Normalized_data/meRIP_RefSeqLongList_Normalized_TKO_2Kb.txt"
TTseq_path <- "/path/4sUDRB/Elongation_rate/Rate_calculation/Elongation_rate_20Kb_Pull_processed_without05.txt"

output <- "/path/4sUDRB/Elongation_rate/Correlations/"
output3 <- "/path/4sUDRB/Elongation_rate/TKOvsWT_rates/"

# -----------
# Open data |
# -----------
m6AWT <- read.table(m6AWT_path, h=T,sep="\t",stringsAsFactors=F)
m6ATKO <- read.table(m6ATKO_path, h=T,sep="\t",stringsAsFactors=F)
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=F,
                    col.names = c("Gene_name","Chr", "Start","End","Exon_number","Strand","TTseq_WT_pull","TTseq_TKO_pull","Size"))
# -------------
# Merge files |
# -------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                list(m6AWT[,c("Chr","Start","End","Strand","meRIP_WT1", "meRIP_WT2","meRIP_WT","Gene_name")],
                     m6ATKO[,c("meRIP_TKO1", "meRIP_TKO2","meRIP_TKO","Gene_name")],
                     TTseq[,c("TTseq_WT_pull", "TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

# ------------------
# Save merged file |
# ------------------
write.table(merged, file = paste0(output,"TTseq_meRIP/TTseq_meRIP.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

# ---------------------------------------------------------------
# Remove outliers from meRIP data and log transform RNApol data |
# ---------------------------------------------------------------
# meRIP_WTmean
Q <- quantile(merged$meRIP_WT, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$meRIP_WT)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_WT <- subset(merged, merged$meRIP_WT > low & merged$meRIP_WT < up)
head(saved_values_WT)
nrow(saved_values_WT)

# For replicate samples use the following:

# meRIP_WT1
#Q <- quantile(merged$meRIP_WT1, probs=c(.25, .75), na.rm = FALSE)
#iqr <- IQR(merged$meRIP_WT1)
#up <-  Q[2]+1.5*iqr # Upper Range  
#low<- Q[1]-1.5*iqr # Lower Range
#saved_values_WT1 <- subset(merged, merged$meRIP_WT1 > low & merged$meRIP_WT1 < up)
#nrow(saved_values_WT1)

# meRIP_WT2
#Q <- quantile(merged$meRIP_WT2, probs=c(.25, .75), na.rm = FALSE)
#iqr <- IQR(merged$meRIP_WT2)
#up <-  Q[2]+1.5*iqr # Upper Range  
#low<- Q[1]-1.5*iqr # Lower Range
#saved_values_WT2 <- subset(merged, merged$meRIP_WT2 > low & merged$meRIP_WT2 < up)
#nrow(saved_values_WT2)

# meRIP_TKOmean
Q <- quantile(merged$meRIP_TKO, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$meRIP_TKO)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_TKO <- subset(merged, merged$meRIP_TKO > low & merged$meRIP_TKO < up)
nrow(saved_values_TKO)

# For replicate samples use the following:

# meRIP_TKO1
#Q <- quantile(merged$meRIP_TKO1, probs=c(.25, .75), na.rm = FALSE)
#iqr <- IQR(merged$meRIP_TKO1)
#up <-  Q[2]+1.5*iqr # Upper Range  
#low<- Q[1]-1.5*iqr # Lower Range
#saved_values_TKO1 <- subset(merged, merged$meRIP_TKO1 > low & merged$meRIP_TKO1 < up)
#nrow(saved_values_TKO1)

# meRIP_TKO2
#Q <- quantile(merged$meRIP_TKO2, probs=c(.25, .75), na.rm = FALSE)
#iqr <- IQR(merged$meRIP_TKO2)
#up <-  Q[2]+1.5*iqr # Upper Range  
#low<- Q[1]-1.5*iqr # Lower Range
#saved_values_TKO2 <- subset(merged, merged$meRIP_TKO2 > low & merged$meRIP_TKO2 < up)
#nrow(saved_values_TKO2)

# ------------------
# Correlation test |
# ------------------
# meRIP WTmean - TTseq WT_pull
png(file = paste0(output, "TTseq_meRIP/Spearman_TTseq_pull_m6A_WTmean.png"))
number <- nrow(saved_values_WT)
ggscatter(saved_values_WT, x = "meRIP_WT", y = "TTseq_WT_pull",
          color="blue",shape = 10,
          size=0.5,
          add.params = list(color = "blue", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = F, cor.method = "spearman",
          ylim=c(0,7.5),
          alpha=.2,
          xlab = "m6A levels", ylab = "Elongation rate (Kb/min)") +
  theme_minimal()+
  stat_cor(method = "spearman", color = "blue", label.x = 0.65, label.y = 6.5,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n m6A levels vs. elongation rates in WT (", number, " genes)"))
dev.off()

# meRIP TKOmean - TTseq TKO_pull
png(file = paste0(output, "TTseq_meRIP/Spearman_TTseq_pull_m6A_TKOmean.png"))
number <- nrow(saved_values_TKO)
ggscatter(saved_values_TKO, x = "meRIP_TKO", y = "TTseq_TKO_pull",
          color="red2",shape = 10,
          size=0.5,
          add.params = list(color = "red2", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = F, cor.method = "spearman",
          ylim=c(0,7.5),
          alpha=.2,
          xlab = "m6A levels", ylab = "Elongation rate (Kb/min)") +
  theme_minimal()+
  stat_cor(method = "spearman", color= "red2", label.x = 0.58, label.y = 6.5,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n m6A levels vs. elongation rates in TKO (", number, " genes)"))
dev.off()

# ----------------------------
# Ratio TKO/WT rates and m6A |
# ----------------------------
TTseq <- merged
TTseq$rates_ratio <- TTseq$TTseq_TKO_pull/TTseq$TTseq_WT_pull
TTseq$m6A_ratio <- TTseq$meRIP_TKO/TTseq$meRIP_WT
TTseq$log_rates <- log(TTseq$rates_ratio)
TTseq <- TTseq[order(TTseq$log_ratio, decreasing = TRUE),]

mean <- mean(TTseq$log_rates)
sd <- sd(TTseq$log_rates)

for(i in 1:nrow(TTseq)) {
  if(TTseq[i,18] >= (mean+2*sd) | TTseq[i,18] <= (mean-2*sd)){
    TTseq[i,19]="??2SD"
  }
  else {
    TTseq[i,19]="95%"
  }
}

colnames(TTseq)[19] <- "CI"

png(file = paste0(output3, "m6Aratio_fast.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_rates>= (mean+2*sd)),]$m6A_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)),]$m6A_ratio)
mean1<-round(mean(TTseq[which(TTseq$log_rates>= (mean+2*sd)),]$m6A_ratio),3)
mean2<-round(mean(TTseq[which(TTseq$log_rates <(mean+2*sd)),]$m6A_ratio),3)
boxplot(TTseq[which(TTseq$log_rates>= (mean+2*sd)),]$m6A_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)),]$m6A_ratio,
        col = c("darkgoldenrod1","darkgoldenrod3"),
        at = c(1,2),
        names = c("TKO fast","Rest"),
        ylab="m6A TKO/WT ratio",
        main = paste0("m6A TKO/WT ratio of TKO faster genes"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2, col= "black")
dev.off()	

png(file = paste0(output3, "m6Aratio_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_rates<= (mean-2*sd)),]$m6A_ratio, TTseq[which(TTseq$log_rates >(mean-2*sd)),]$m6A_ratio)
mean1<-round(mean(TTseq[which(TTseq$log_rates<= (mean-2*sd)),]$m6A_ratio),3)
mean2<-round(mean(TTseq[which(TTseq$log_rates >(mean-2*sd)),]$m6A_ratio),3)
boxplot(TTseq[which(TTseq$log_rates<= (mean-2*sd)),]$m6A_ratio, TTseq[which(TTseq$log_rates >(mean-2*sd)),]$m6A_ratio,
        col = c("lightcyan3","lightcyan4"),
        at = c(1,2),
        names = c("TKO slow","Rest"),
        ylab="m6A TKO/WT ratio",
        main = paste0("m6A TKO/WT ratio of TKO slower genes"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2+.01, col= "black")
dev.off()	

png(file = paste0(output3, "m6Aratio_fast_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_rates>= (mean+2*sd) | TTseq$log_rates<= (mean-2*sd)),]$m6A_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)|TTseq$log_rates >(mean-2*sd)),]$m6A_ratio)
mean1<-round(mean(TTseq[which(TTseq$log_rates>= (mean+2*sd) | TTseq$log_rates<= (mean-2*sd)),]$m6A_ratio),3)
mean2<-round(mean(TTseq[which(TTseq$log_rates <(mean+2*sd)|TTseq$log_rates >(mean-2*sd)),]$m6A_ratio),3)
boxplot(TTseq[which(TTseq$log_rates>= (mean+2*sd) | TTseq$log_rates<= (mean-2*sd)),]$m6A_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)|TTseq$log_rates >(mean-2*sd)),]$m6A_ratio,
        col = c("springgreen3","springgreen4"),
        at = c(1,2),
        names = c("TKO fast & slow","Rest"),
        ylab="m6A TKO/WT ratio",
        main = paste0("m6A TKO/WT ratio of TKO faster and slower genes"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2+.01, col= "black")
dev.off()	
