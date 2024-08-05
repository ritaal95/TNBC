#install.packages("corrplot")
#install.packages("readxl")
#install.packages("tidyverse")

#Select the Directory where you can find the RData with the information regarding the TCGA-TNBC cohort
#obtained from the "Select the TNBC.R" code

#############################################################################################################
######################################### Immune cells frequency: ###########################################
######################### Comparing between high and low ST6GALNAC1 expression groups #######################
###################################### Correlations with ST6GALNAC1 #########################################
#############################################################################################################

load("data/logTPMs.RData")

library(readxl)
library(corrplot)
library(tidyverse)

####Separating samples with the median to define high and low####
Exp_level<-function(Gene){
  Info <- list()
  for (i in 1:length(Gene)){
    GeneExpLevels <- logTPMs[Gene[i],]
    GeneExpStatus <- ifelse(GeneExpLevels > median(GeneExpLevels), "High", "Low")
    Final_status <- list(data.frame(names(GeneExpLevels),GeneExpLevels,GeneExpStatus))
    Info <- append(Info,Final_status)
  }
  names(Info) <- Gene
  return (Info)
}
#Name the gene of interest
Gene <- c("ST6GALNAC1")
GenesNames <- c("CD44", "CDH2", "CTNNB1", "MMP9", "POU5F1", "SOX2", "MKI67","MYC", "SALL4")

GeneExp <- c()
for (i in GenesNames){
  GeneExp[[i]] <- logTPMs[i,]
}
GeneExp<-data.frame(GeneExp)
GeneExp <- GeneExp[sort(rownames(GeneExp)),]
Info <- as.data.frame(Exp_level(Gene))
Info <- Info[,-1]
colnames(Info) <- c("ST6GALNAC1", "Group")
Info <- Info[sort(rownames(Info)),]
GeneTable <- cbind(Info, GeneExp)
GeneTable <- GeneTable[,c(2,1,3,4,5,6,7,8,9,10,11)]

##Acquiring the table with the cell type fractions in each sample
data <- as.data.frame(read_excel('data/Filtrates.xlsx'))
rownames(data) <- data$Patients
data <- data[,-1]
colnames(data)
colnames(data) <- c("Uncharacterized Cells", "Treg Cells", "Neutrophiles", "NK Cells", "Monocytes", "Macrophages M1", "Macrophages M2", "Dendritic Cells", "CD8 T Cells", "CD4 T Cells", "B Cells")
data <- data[order(row.names(data)), ]
GeneTable <- GeneTable[order(row.names(GeneTable)), ]
data <- cbind(GeneTable$ST6GALNAC1, data)
names(data)[names(data) == 'GeneTable$ST6GALNAC1'] <- 'ST6GALNAC1'

###Without Uncharacterized cells
options(scipen = 999, digits = 3)
Test <- colnames(data[,3:12])
Cor <- matrix(data=0, nrow=length(Test), ncol=2)
rownames(Cor) <- Test
colnames(Cor) <- c("CorValue", "pValue")
for(i in 1:length(Test)){
  cur_result <- cor.test(data[,"ST6GALNAC1"], data[,Test[i]], method="spearman", exact=F)
  Cor[i,] <- c(cur_result$estimate, cur_result$p.value)
}

adj.pValues <- p.adjust(Cor[,"pValue"], method="fdr")
Cor <- cbind(Cor, adj.pValues)
Cor<- as.data.frame(Cor)
Cor$CorValue <- as.numeric(as.character(Cor$CorValue))
Cor$pValue <- as.numeric(as.character(Cor$pValue))
Cor$adj.pValues <- as.numeric(as.character(Cor$adj.pValues))
TCGACor <- Cor
write.table(TCGACor, file = "TCGACor.csv")

M <- as.matrix(Cor$CorValue)
colnames(M) <- "ST6GALNAC1"
rownames(M) <- rownames(Cor)
p.mat <- as.matrix(Cor$adj.pValues)
colnames(p.mat) <- "ST6GALNAC1"
rownames(p.mat) <- rownames(Cor)
#horizontal heatmap
M <- t(M)
p.mat <- t(p.mat)

png("HeatmapSI.png", res=600, width=10100, height=3500)
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
         mar = c(0, 5, 0, 14)
)
colorlegend(
  colorRampPalette(c("#2166AC","white","#D6604D"))(10),
  c(seq(-1,1,0.2)),
  xlim=c(0.5,10.5),
  ylim=c(-0.2,0.4),
  align="l",
  vertical=FALSE,
  addlabels=TRUE,
  cex=2
)
text(2.2, -0.5, "Spearman Correlation", col = "black", cex = 2)
dev.off()

GeneExpLevels <- data$ST6GALNAC1
GeneExpStatus <- ifelse(GeneExpLevels > median(GeneExpLevels), "High", "Low")
Group <- GeneExpStatus
data <- cbind(Group, data)

high_subset <- data[data$Group=='High', 4:13]
low_subset <- data[data$Group=='Low', 4:13]
high_subset <- high_subset %>% mutate_if(is.character, as.numeric)
low_subset <- low_subset %>% mutate_if(is.character, as.numeric)

cellType <- colnames(data[,4:13])

Cor <- matrix(data=0, nrow=length(cellType), ncol=1)
rownames(Cor) <- cellType
colnames(Cor) <- c("pValue")
for (i in cellType) {
  wilcoxTest <- wilcox.test(high_subset[,i],low_subset[,i], paired = FALSE )
  Cor[i,] <- c(wilcoxTest$p.value)
}
adj.pValues <- p.adjust(Cor[,"pValue"], method="fdr")
Cor <- cbind(Cor, adj.pValues)

#Cell type pValue, adj. pValue
#Macrophages M2 0.00000316, 0.0000316
#Treg Cells 0.00003426, 0.0001713
#B Cells 0.00008090, 0.0002697
#CD4 T Cells 0.00138362, 0.0034590
#CD8 T Cells 0.03303383, 0.0660677
#Monocytes 0.08671830, 0.1445305
#Macrophages M1 0.21542079, 0.3077440
#NK Cells 0.58737381, 0.7342173
#Neutrophiles 0.75744004, 0.8131563
#Dendritic Cells 0.81315635, 0.8131563
Cor <- as.data.frame(Cor)
write.table(Cor, file = "WilcoxTes_ST6GALNAC1_immune.csv")

#prepare data for graphpad plot
DataLow <- matrix(data=0, nrow=length(cellType), ncol=3)
rownames(DataLow) <- cellType
colnames(DataLow) <- c("Mean", "SD", "N")
for(i in cellType){
  meanLow <- mean(low_subset[,i], na.rm=T)
  sdLow <- sd(low_subset[,i], na.rm=T)
  nLow<- length(low_subset[,i])
  DataLow[i,] <- c(meanLow, sdLow, nLow)
}

DataHigh <- matrix(data=0, nrow=length(cellType), ncol=3)
rownames(DataHigh) <- cellType
colnames(DataHigh) <- c("Mean", "SD", "N")
for(i in cellType){
  meanHigh <- mean(high_subset[,i], na.rm=T)
  sdHigh <- sd(high_subset[,i], na.rm=T)
  nHigh<- length(high_subset[,i])
  DataHigh[i,] <- c(meanHigh, sdHigh, nHigh)
}

write.table(DataHigh, "High stats ST6GALNAC1.csv", sep=";", quote=F, col.names=NA)
write.table(DataLow, "Low stats ST6GALNAC1.csv", sep=";", quote=F, col.names=NA)

############################################################## END ##########################################
