#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("EnsDb.Hsapiens.v86")
#install.packages(corrplot)
#install.packages("openxlsx")

#Select the Directory where you can find the RData with the information regarding the TCGA-TNBC cohort
#obtained from the "Select the TNBC.R" code

#############################################################################################################
######################### Correlations between ST6GALNAC1 and the TGF-B pathway genes #######################
#############################################################################################################

#Load necessary libraries
library("EnsDb.Hsapiens.v86")
library(corrplot)
library(openxlsx)

#Load Data
load("data/TNBC_data.RData")

#Normality of data
options(scipen = 999, digits = 3)
geneIDs <- c("SMAD2", "SMAD3", "SMAD4", "SMAD7", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "SMURF1", "SMURF2", "ZFYVE9", "ZFYVE16")
myGenes <- geneIDs
normality <- matrix(data=0, nrow=length(myGenes), ncol=2)
rownames(normality) <- myGenes
colnames(normality) <- c("StatisticW", "pValue")
for(i in 1:length(myGenes)){
  normality_result <- shapiro.test(logTPMs[myGenes[i],])
  normality [i,] <- c(normality_result$statistic, normality_result$p.value)
}
#not all are normal so I'll use the non-parametric spearman correlation test

#Correlation between ST6GALNAC1 and the biomarkers genes: 
options(scipen = 999, digits = 3)
myGenes <- geneIDs
Cor <- matrix(data=0, nrow=length(myGenes), ncol=2)
rownames(Cor) <- myGenes
colnames(Cor) <- c("CorValue", "pValue")
for(i in 1:length(myGenes)){
  cur_result <- cor.test(logTPMs["ST6GALNAC1",], logTPMs[myGenes[i],], method="spearman", exact=F)
  Cor[i,] <- c(cur_result$estimate, cur_result$p.value)
}

adj.pValues <- p.adjust(Cor[,"pValue"], method="fdr")
Cor <- cbind(Cor, adj.pValues)
geneNames <- genes(EnsDb.Hsapiens.v86, return.type="data.frame"); 
Cor<- as.data.frame(Cor)
Cor$CorValue <- as.numeric(as.character(Cor$CorValue))
Cor$pValue <- as.numeric(as.character(Cor$pValue))
Cor$adj.pValues <- as.numeric(as.character(Cor$adj.pValues))
index_SignCorrelations <- which(Cor[,"adj.pValues"] < 0.05)
length(index_SignCorrelations)
signCorrelations_AllGenes <- Cor[index_SignCorrelations,]

Cor <- Cor[,-c(2)]
TCGACor <- Cor
write.xlsx(TCGACor, file = "TCGACor.xlsx", rownames = T)

#Prepare data to obtain the heatmap
M <- as.matrix(TCGACor$CorValue)
colnames(M) <- "ST6GALNAC1"
rownames(M) <- rownames(TCGACor)
p.mat <- as.matrix(TCGACor$adj.pValues)
colnames(p.mat) <- "ST6GALNAC1"
rownames(p.mat) <- rownames(TCGACor)

png("Heatmap ST6GALNAC1_TGFB.png",res=300, width=2000, height=2700)
corrplot(M, 
         method="color", 
         type="full", 
         cl.pos = "n",
         tl.srt=65,
         order="original", 
         p.mat = p.mat, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 1.5,
         diag=TRUE, 
         tl.col="black", 
         tl.cex =1.75,
         col=colorRampPalette(c("#2166AC","white","#D6604D"))(200),
         addgrid.col="black"
)
colorlegend(
  colorRampPalette(c("#2166AC","white","#D6604D"))(20),
  c(seq(-1,1,0.25)),
  xlim=c(1.65,2.25),
  ylim=c(0.5,14.5),
  align="l",
  vertical=TRUE,
  addlabels=TRUE,
  cex=1
)
text(3.3, 3, "Spearman Correlation", col = "black", cex = 1.5, srt=90)

dev.off()


