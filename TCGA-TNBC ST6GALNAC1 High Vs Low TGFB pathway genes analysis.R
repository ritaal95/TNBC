#############################################################################################################
########### Comparing TGFB pathway gene expression between high and low ST6GALNAC1 expression groups#########
#############################################################################################################

##Libraries used
#install.packages(dplyr)

##Load libraries
#library(dplyr)

##Select the Directory where you can find the RData with the information regarding the TCGA-TNBC cohort obtained from the "Select the TNBC.R" code

##Load the gene expression data
load("data/logTPMs.RData")

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
GenesNames <- c("SMAD2", "SMAD3", "SMAD4", "SMAD7", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "SMURF1", "SMURF2", "ZFYVE9", "ZFYVE16")

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

high_subset <- GeneTable[GeneTable$Group=='High', 3:16]
low_subset <- GeneTable[GeneTable$Group=='Low', 3:16]
high_subset <- high_subset %>% mutate_if(is.character, as.numeric)
low_subset <- low_subset %>% mutate_if(is.character, as.numeric)


Cor <- matrix(data=0, nrow=length(GenesNames), ncol=1)
rownames(Cor) <- GenesNames
colnames(Cor) <- c("pValue")
for (i in GenesNames) {
  wilcoxTest <- wilcox.test(high_subset[,i],low_subset[,i], paired = FALSE )
  Cor[i,] <- c(wilcoxTest$p.value)
}
adj.pValues <- p.adjust(Cor[,"pValue"], method="fdr")
Cor <- cbind(Cor, adj.pValues)
#Gene pValue, adj pValue
#POU5F1 9.261877e-05, 0.0004167845
#MYC 4.710460e-05, 0.0004167845
#SOX2 3.489783e-02, 0.1046935019
#MKI67 1.318950e-01, 0.2967638156
#CDH2 2.494203e-01, 0.3762386271 
#SALL4 2.508258e-01, 0.3762386271
#CTNNB1 6.291771e-01, 0.8089419331
#CD44 9.959157e-01, 0.9959157048
#MMP9 9.850251e-01, 0.9959157048

##Obtain data do make the plots in graphPad prims
DataLow <- matrix(data=0, nrow=length(GenesNames), ncol=3)
rownames(DataLow) <- GenesNames
colnames(DataLow) <- c("Mean", "SD", "N")
for(i in 1:length(GenesNames)){
  meanLow <-mean(as.numeric(low_subset[,GenesNames[i]], na.rm=T))
  sdLow <- sd(as.numeric(low_subset[,GenesNames[i]], na.rm=T))
  nLow <- length(low_subset[,GenesNames[i]])
  DataLow[i,] <- c(meanLow, sdLow, nLow)
}

DataHigh <- matrix(data=0, nrow=length(GenesNames), ncol=3)
rownames(DataHigh) <- GenesNames
colnames(DataHigh) <- c("Mean", "SD", "N")
for(i in 1:length(GenesNames)){
  meanHigh <-mean(as.numeric(high_subset[,GenesNames[i]], na.rm=T))
  sdHigh <- sd(as.numeric(high_subset[,GenesNames[i]], na.rm=T))
  nHigh <- length(high_subset[,GenesNames[i]])
  DataHigh[i,] <- c(meanHigh, sdHigh, nHigh)
}
Cor <- as.data.frame(Cor)

#Save the data tables
write.table(DataHigh, "S6GALNAC1 High stats TGFB.csv", sep=";", quote=F, col.names=NA)
write.table(DataLow, "S6GALNAC1 Low stats TGFB.csv", sep=";", quote=F, col.names=NA)
write.table(Cor, "ST6GALNAC1 groups.csv")
write.table(high_subset, "ST6GALNAC1 High Group TGFB.csv")
write.table(low_subset, "ST6GALNAC1 Low Group TGFB.csv")

############################################################## END ##########################################
