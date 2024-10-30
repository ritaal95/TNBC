#############################################################################################################
########################### Comparing biomarkers expression between STn- and STn+  groups####################
#############################################################################################################

##Librarie used
#install.packages(openxlsx)

##Load librarie
#library(openxlsx)

##Select the Directory where you can find the Table S2 with the information regarding the HSJ-TNBC cohort
TNBC_Data <- as.data.frame(read.xlsx("data/HSJ-TNBC_Data.xlsx"))
Biomarkers <- TNBC_Data[,24:37]
rownames(Biomarkers) <- TNBC_Data$ID

X <- split(Biomarkers, Biomarkers$`Sialyl-Tn_CM_Tumor`)
Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:14]) 
STnNegative <- Y[[1]]
STnPositive <- rbind(Y[[2]], Y[[3]], Y[[4]], Y[[5]], Y[[6]], Y[[7]], Y[[8]])

STnPositive$CD44_CM_Tumor <- as.numeric(STnPositive[,1])
STnPositive$CD44_CM_Stroma <- as.numeric(STnPositive[,2])
STnPositive$`N-Cad_CM_tumor` <- as.numeric(STnPositive[,3])
STnPositive$`N-Cad_N_Tumor` <- as.numeric(STnPositive[,4])
STnPositive$`B-Cat_CM_Tumor` <- as.numeric(STnPositive[,5])
STnPositive$MMP9_C_Tumor <- as.numeric(STnPositive[,6])
STnPositive$MMP9_C_Stroma <- as.numeric(STnPositive[,7])
STnPositive$Oct4_N_Tumor <- as.numeric(STnPositive[,8])
STnPositive$Oct4_C_Tumor <- as.numeric(STnPositive[,9])
STnPositive$Sox2_N_Tumor <- as.numeric(STnPositive[,10])
STnPositive$KI67_N_Tumor <- as.numeric(STnPositive[,11])
STnPositive$MYC_N_Tumor <- as.numeric(STnPositive[,12])
STnPositive$Sall4_N_Tumor <- as.numeric(STnPositive[,13])
STnPositive$`Sialyl-Tn_CM_Tumor` <- as.numeric(STnPositive[,14])

STnNegative$CD44_CM_Tumor <- as.numeric(STnNegative[,1])
STnNegative$CD44_CM_Stroma <- as.numeric(STnNegative[,2])
STnNegative$`N-Cad_CM_tumor` <- as.numeric(STnNegative[,3])
STnNegative$`N-Cad_N_Tumor` <- as.numeric(STnNegative[,4])
STnNegative$`B-Cat_CM_Tumor` <- as.numeric(STnNegative[,5])
STnNegative$MMP9_C_Tumor <- as.numeric(STnNegative[,6])
STnNegative$MMP9_C_Stroma <- as.numeric(STnNegative[,7])
STnNegative$Oct4_N_Tumor <- as.numeric(STnNegative[,8])
STnNegative$Oct4_C_Tumor <- as.numeric(STnNegative[,9])
STnNegative$Sox2_N_Tumor <- as.numeric(STnNegative[,10])
STnNegative$KI67_N_Tumor <- as.numeric(STnNegative[,11])
STnNegative$MYC_N_Tumor <- as.numeric(STnNegative[,12])
STnNegative$Sall4_N_Tumor <- as.numeric(STnNegative[,13])
STnNegative$`Sialyl-Tn_CM_Tumor` <- as.numeric(STnNegative[,14])


Test <- colnames(Biomarkers)[1:13]
Cor <- matrix(data=0, nrow=length(Test), ncol=1)
rownames(Cor) <- Test
colnames(Cor) <- c("pValue")
for(i in 1:length(Test)){
  wilcoxTest <- wilcox.test(STnNegative[,Test[i]], STnPositive[,Test[i]])
  Cor[i,] <- c(wilcoxTest$p.value)
}
adj.pValues <- p.adjust(Cor[,"pValue"], method="fdr")
Cor <- cbind(Cor, adj.pValues)
Cor <- as.data.frame(Cor)
#Biomarker pValue, adj pValue
#MYC_N_Tumor 0.002044194, 0.02657453
#KI67_N_Tumor 0.018820096, 0.12233063
#MMP9_C_Tumor 0.050270745, 0.21783989
#CD44_CM_Stroma 0.302682128, 0.49348355
#N-Cad_CM_tumor 0.303682184, 0.49348355
#MMP9_C_Stroma 0.159326913, 0.49348355
#Oct4_N_Tumor 0.279905979, 0.49348355
#Oct4_C_Tumor 0.302271338, 0.49348355
#B-Cat_CM_Tumor 0.419047474, 0.54476172 
#Sox2_N_Tumor 0.412487825, 0.54476172
#CD44_CM_Tumor 0.574495095, 0.67894875
#N-Cad_N_Tumor 0.697983222, 0.75614849
#Sall4_N_Tumor 0.997259081, 0.99725908
write.xlsx(Cor, "WilcoxTes_STn+VsSTn-.xlsx")

##Obtain data do make the plots in graphPad prims
Test <- colnames(Biomarkers)[1:13]
DataNeg <- matrix(data=0, nrow=length(Test), ncol=3)
rownames(DataNeg) <- Test
colnames(DataNeg) <- c("Mean", "SD", "N")
for(i in 1:length(Test)){
  meanSTnNeg <- mean(STnNegative[,Test[i]], na.rm=T)
  sdSTnNeg <- sd(STnNegative[,Test[i]], na.rm=T)
  nNeg <- length(STnNegative[,Test[i]])
  DataNeg[i,] <- c(meanSTnNeg, sdSTnNeg, nNeg)
}

Test <- colnames(Biomarkers)[1:13]
DataPos <- matrix(data=0, nrow=length(Test), ncol=3)
rownames(DataPos) <- Test
colnames(DataPos) <- c("Mean", "SD", "N")
for(i in 1:length(Test)){
  meanSTnPos <- mean(STnPositive[,Test[i]], na.rm=T)
  sdSTnPos <- sd(STnPositive[,Test[i]], na.rm=T)
  nPos <- length(STnPositive[,Test[i]])
  DataPos[i,] <- c(meanSTnPos, sdSTnPos, nPos)
}

#Save the data tables
write.table(DataPos, "STn+ stats.csv", sep=";", quote=F, col.names=NA)
write.table(DataNeg, "STn- stats.csv", sep=";", quote=F, col.names=NA)

