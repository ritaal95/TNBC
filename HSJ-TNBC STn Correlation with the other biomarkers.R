#install.packages(corrplot)
#install.packages(openxlsx)

#Select the Directory where you can find the table with the information regarding the HSJ-TNBC cohort
#Table S1 from Supplementary File S2

#############################################################################################################################
##################################### STn correlation with the other biomarkers analysed ####################################
#############################################################################################################################

library(openxlsx)
library(corrplot)

TNBC_Data <- as.data.frame(read.xlsx("data/HSJ-TNBC_Data.xlsx"))

#Create subset tables with the clinical informationa and the biomarkers information
ClinicalInfo <- TNBC_Data[,2:23]
rownames(ClinicalInfo) <- TNBC_Data$ID
Biomarkers <- TNBC_Data[,24:37]
rownames(Biomarkers) <- TNBC_Data$ID


#Testing for normality
Test <- colnames(Biomarkers)[1:14]
normality <- matrix(data=0, nrow=length(Test), ncol=2)
rownames(normality) <- Test
colnames(normality) <- c("StatisticW", "pValue")
for(i in 1:length(Test)){
  normality_result <- shapiro.test(Biomarkers[,i])
  normality [i,] <- c(normality_result$statistic, normality_result$p.value)
}
#data is not normal so apply the nonparametric Spearman corelation

#Correlation of all markers with the STn
Test <- colnames(Biomarkers)[1:13]
Cor <- matrix(data=0, nrow=length(Test), ncol=2)
rownames(Cor) <- Test
colnames(Cor) <- c("CorValue", "pValue")
for(i in 1:length(Test)){
  cur_result <- cor.test(as.numeric(Biomarkers[,14]), as.numeric(Biomarkers[,Test[i]]), method="spearman", exact = FALSE)
  Cor[i,] <- c(cur_result$estimate, cur_result$p.value)
}
adj.pValues <- p.adjust(Cor[,"pValue"], method="fdr")
Cor <- cbind(Cor, adj.pValues)
Cor<- as.data.frame(Cor)
Cor$CorValue <- as.numeric(as.character(Cor$CorValue))
Cor$pValue <- as.numeric(as.character(Cor$pValue))
Cor$adj.pValues <- as.numeric(as.character(Cor$adj.pValues))
CohortCor <- Cor
write.xlsx(CohortCor, file = "HSJCor.xlsx", rownames=T)
#prepare data to make a hetmap
rownames(CohortCor) <- c("CD44 CMT", "CD44 CMS", "N-Cad CMT", "N-Cad NT", "β-Cat CMT", "MMP9 CT", "MMP9 CS", "Oct4 NT", "Oct4 CT", "Sox2 NT", "Ki67 NT", "c-Myc NT", "Sall4 NT")

#vertical plot
M <- as.matrix(CohortCor$CorValue)
colnames(M) <- "STn CMT"
rownames(M) <- rownames(CohortCor)
p.mat <- as.matrix(CohortCor$adj.pValues)
colnames(p.mat) <- "STn CMT"
rownames(p.mat) <- rownames(CohortCor)
M <- t(M)
p.mat <- t(p.mat)

png("Heatmap STn.png", res=600, width=8000, height=3000)
corrplot(M, 
         method="color", 
         type="full", 
         cl.pos = "n",
         tl.srt=65,
         order="original", 
         p.mat = p.mat, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 1.75,
         diag=TRUE, 
         tl.col="black", 
         tl.cex =1.75,
         col=colorRampPalette(c("#2166AC","white","#D6604D"))(20),
         addgrid.col="black"
)
colorlegend(
  colorRampPalette(c("#2166AC","white","#D6604D"))(20),
  c(seq(-1,1,0.25)),
  xlim=c(0.5,13.5),
  ylim=c(-0.2,0.4),
  align="l",
  vertical=FALSE,
  addlabels=TRUE,
  cex=1.75
)
text(1.8, -0.5, "Spearman Correlation", col = "black", cex = 1.5)

dev.off()

############################################################## END ##########################################################