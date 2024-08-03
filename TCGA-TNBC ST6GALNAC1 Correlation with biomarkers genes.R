#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("EnsDb.Hsapiens.v86")
#install.packages(corrplot)
#install.packages(openxlsx)

#Select the Directory where you can find the RData with the information regarding the TCGA-TNBC cohort
#obtained from the "Select the TNBC.R" code

#############################################################################################################
####################### Correlations between ST6GALNAC1 and the other biomarkers genes ######################
#############################################################################################################

#Load necessary libraries
library("EnsDb.Hsapiens.v86")
library(corrplot)
library(openxlsx)

#Load Data
load("data/TNBC_data.RData")
load("data/logTPMs.RData")

#Normality of data
options(scipen = 999, digits = 3)
geneIDs <- c("CD44","CDH2", "CTNNB1", "MMP9","POU5F1", "SOX2", "MKI67",  "MYC", "SALL4")
myGenes <- geneIDs
normality <- matrix(data=0, nrow=length(myGenes), ncol=2)
rownames(normality) <- myGenes
colnames(normality) <- c("StatisticW", "pValue")
for(i in 1:length(myGenes)){
  normality_result <- shapiro.test(logTPMs[myGenes[i],])
  normality [i,] <- c(normality_result$statistic, normality_result$p.value)
}
#not all data is normal so apply the nonparametric Spearman corelation

#Correlation between ST6GALNAC1 and the biomarkers genes: 
options(scipen = 999, digits = 3)
geneIDs <- c("CD44","CDH2", "CTNNB1", "MMP9","POU5F1", "SOX2", "MKI67",  "MYC", "SALL4")
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
TCGACor <- Cor
write.xlsx(TCGACor, file = "TCGACor.xlsx", rownames = T)

#Prepare data to obtain the heatmap
M <- as.matrix(TCGACor$CorValue)
colnames(M) <- "ST6GALNAC1"
rownames(M) <- rownames(TCGACor)
p.mat <- as.matrix(TCGACor$adj.pValues)
colnames(p.mat) <- "ST6GALNAC1"
rownames(p.mat) <- rownames(TCGACor)

#horizontal heatmap
M <- t(M)
p.mat <- t(p.mat)

png("Heatmap ST6GALNAC1.png", res=600, width=10100, height=3500)
corrplot(M, 
         method="color", 
         type="full", 
         cl.pos = "n",
         tl.srt=65,
         order="original", 
         p.mat = p.mat, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 2,
         diag=TRUE, 
         tl.col="black", 
         tl.cex =2,
         col=colorRampPalette(c("#2166AC","white","#D6604D"))(10),
         addgrid.col="black",
         mar = c(0, 7, 0, 16.85)
)
colorlegend(
  colorRampPalette(c("#2166AC","white","#D6604D"))(10),
  c(seq(-1,1,0.2)),
  xlim=c(0.5,9.5),
  ylim=c(-0.2,0.4),
  align="l",
  vertical=FALSE,
  addlabels=TRUE,
  cex=2
)
text(2, -0.5, "Spearman Correlation", col = "black", cex = 2)

dev.off()


############################################################## END ##########################################
