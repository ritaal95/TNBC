#############################################################################################
################ Contingency tables to analyse STn with the clinical fetures ################
############################### Supplementary Material Nr XX ################################
#############################################################################################

##Libraries Used
#install.packages(tidyverse)
#install.packages(readxl)

##Load libraries
#library(readxl)
#library(tidyverse)

##Select the Directory where you can find the Table S2 with the information regarding the HSJ-TNBC cohort
##Preparing the data table for analysis
TNBC_Data <- as.data.frame(read_xlsx("data/HSJ-TNBC_Data.xlsx"))
ClinicalInfo <-  dplyr::select(TNBC_Data, c("Laterality","Age","Size",
                                            "Grade","Margin","Invasion",
                                            "SLN_Met","GG_Met","Dist_Met_Dt",
                                            "Relapse","Sialyl-Tn_CM_Tumor"))
rownames(ClinicalInfo) <- TNBC_Data$ID

#Separating the samples in negative and positive for STn
STnExpLevels<-ClinicalInfo$`Sialyl-Tn_CM_Tumor`
STnExpStatus<-ifelse(STnExpLevels < "1", "STn-", "STn+")

##Testing Age versus STn expression
#Vector with patient's age
Age<-as.vector(ClinicalInfo$Age)
#Separating the samples in groups by the median age
Age <- ifelse(Age >= median(Age), "≥ 62y", "< 62y")
#creating a contingency table with the age groups and the STn groups 
#statistical analysis using the Fisher exact test
contTable<-table(STnExpStatus, Age)
fisher.test(contTable)
##p-value = 0.02046
##odds ratio 2.97382
#ploting the contingency table and saving it as an image
plotdata<-apply(contTable,1, function(x)x/sum(x)*100)
png("Age_contTable.png", res = 1200, width = 5000, height = 4000)
#grDevices::windows()
par(xpd=F, mar=par()$mar+c(-2,0,-1,3.8))
barplot(plotdata, 
        legend.text = T, 
        col = c("grey3", "grey83"), 
        args.legend = list(x = "topright",inset = c(- 0.4, 0))
)
mtext("p = 0.020", line = 0.5)
dev.off()

##Testing tumour Size versus STn expression
#Vector with patient's tumour size
Size<-as.vector(ClinicalInfo$Size)
#Separating the samples in groups by the median size of the tumours
Size <- ifelse(Size<25, "< 25mm", "≥ 25mm")
#creating a contingency table with the size groups and the STn groups
#statistical analysis using the Fisher exact test
contTable<-table(STnExpStatus, Size)
fisher.test(contTable)
##p-value = 0.09093
##odds ratio 2.217535
#ploting the contingency table and saving it as an image
plotdata<-apply(contTable,1, function(x)x/sum(x)*100)
png("Size_contTable.png", res = 1200, width = 5000, height = 4000)
par(xpd=F, mar=par()$mar+c(-2,0,-1,3.8))
barplot(plotdata, 
        legend.text = T, 
        col = c("grey3", "grey83"), 
        args.legend = list(x = "topright",inset = c(- 0.48, 0))
)
mtext("p = 0.090", line=0.5)
dev.off()

##Testing tumour Grade versus STn expression
#Vector with patient's tumour grade
Grade<-as.vector(ClinicalInfo$Grade)
#Separating the samples in groups by the Grade in "Grade I/II" and "Grade III"
tumorGrade<-ifelse(Grade=="0"|Grade=="1", "Grade I/II", "Grade III")
#creating a contingency table with the Grade groups and the STn groups
#statistical analysis using the Fisher exact test
contTable<-table(STnExpStatus, tumorGrade)
fisher.test(contTable)
##p-value = 0.01201
##odds ratio 0.2799213
#ploting the contingency table and saving it as an image
plotdata<-apply(contTable,1, function(x)x/sum(x)*100)
png("Grade_contTable.png", res = 1200, width = 5000, height = 4000)
par(xpd=F, mar=par()$mar+c(-2,0,-1,3.8))
barplot(plotdata,
        legend.text = T, 
        col = c("grey3", "grey83"), 
        args.legend = list(x = "topright",inset = c(- 0.53, 0))
)
mtext("p = 0.012", line = 0.5)
dev.off()

########################################### END #############################################

