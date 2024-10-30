#############################################################################################################################
######################################### Survival Analysis of HSJ-TNBC Patients Samples ########################################
#############################################################################################################################

##Libraries used
#install.packages(survival)
#install.packages(survminer)
#install.packages(readxl)

##Load libraries
#library(survival)
#library(survminer)
#library(readxl)

##Select the Directory where you can find the Table S2 with the information regarding the HSJ-TNBC cohort
##Preparing the data table for analysis
TNBC_Data <- as.data.frame(read_xlsx("data/HSJ-TNBC_Data.xlsx"))
ClinicalInfo <-  dplyr::select(TNBC_Data, c("Dx_Date", "Dist_Met_Dt", "Relapse", 
                                            "Relapse_Dt", "Alive", "Death", 
                                            "Age", "Size", 
                                            "Grade", "SLN_Met", "GG_Met", "Invasion"))
rownames(ClinicalInfo) <- TNBC_Data$ID
Last_followup <- ifelse(ClinicalInfo$Alive == 1, "2015-01-01", NA)
ClinicalInfo <- cbind(ClinicalInfo, Last_followup)
days_to_last_followup <- as.numeric(difftime(ClinicalInfo$Last_followup, ClinicalInfo$Dx_Date))
days_to_death <- as.numeric(difftime(ClinicalInfo$Death, ClinicalInfo$Dx_Date))
ClinicalInfo <- cbind(ClinicalInfo, days_to_death)
ClinicalInfo <- cbind(ClinicalInfo, days_to_last_followup)
death_status <- ifelse(ClinicalInfo$Alive == 1, "0", "1")
ClinicalInfo <- cbind(ClinicalInfo, death_status)
ClinicalInfo$Dx_Date <- as.Date(ClinicalInfo$Dx_Date)
ClinicalInfo$Dist_Met_Dt <- as.Date(ClinicalInfo$Dist_Met_Dt)
ClinicalInfo$Relapse_Dt <- as.Date(ClinicalInfo$Relapse_Dt)
ClinicalInfo$Last_followup <- as.Date(ClinicalInfo$Last_followup)
ClinicalInfo$Death <- as.Date(ClinicalInfo$Death)

##Obtaining the Progression Free Survival of our TNBC Cohort
#Creating the Progression Free Survival (PFS) and the event columns to obtain the Progression free survival 
PFS <- as.data.frame(ifelse(is.na(pmin(ClinicalInfo[,2],ClinicalInfo[,4],ClinicalInfo[,6], na.rm = TRUE)-ClinicalInfo[,1]),
                            ClinicalInfo[,15], 
                            pmin(ClinicalInfo[,2],ClinicalInfo[,4], ClinicalInfo[,6],na.rm = TRUE)-ClinicalInfo[,1]))
colnames(PFS) <- "PFS"
event <- as.data.frame(ifelse(is.na(ClinicalInfo[,2]), 
                              ifelse(is.na(ClinicalInfo[,4]), 
                                     ifelse(is.na(ClinicalInfo[,6]), 
                                            0, 1), 1), 1))
colnames(event) <- "event"
ClinicalInfo <- cbind(ClinicalInfo, PFS, event)

#creating a table named data with the time in years and with the status that in this case is the occurrence of an event
#event -> recurrence and or death as 1 if it happen and 0 if it didn't happen
time <- apply(ClinicalInfo[c("PFS")], 1, function(x) as.numeric(na.omit(x)))/365
Status <- as.numeric(unlist(ClinicalInfo[,"event"]))
data <- data.frame("time"=time, "status"=Status)
colnames(data)<-c('time','status')

#creating a survival object and the survival curve
surv_data <- Surv(data$time, data$status)
surv_fit <- survfit(surv_data~1, data=data)
summary(surv_fit)
survfit(formula = surv_data ~ 1, data = data)
#50% probability of progression free survival reach after 8.5 years (8 years)
min(data$time)
#0
max(data$time)
#14.03, 14 years

#plotting the survival curve and saving it as an image
png("PFS.png", res=1200, width=7000, height=7000)
ggsurv<-ggsurvplot(surv_fit,
           xlab = "Time (years)",
           ylab = "Progression-Free Survival Probability",
           palette = c("red"),
           conf.int = T,
           pval.coord = c(0,0),
           risk.table = FALSE,
           risk.table.col = " ",
           risk.table.height = 0.25,
           ggtheme = theme_bw(base_family =  "sans"),
           censor = T,
           data=data,
           font.title=c(20, "bold", "black"),
           font.x=c(15, "bold", "black"), 
           font.y=c(15, "bold", "black"), 
           font.tickslab=c(15, "bold", "black"),
           font.legend=c(13, "bold", "black"),
           surv.median.line = "hv"
)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 7.6, y = 0.47, # x and y coordinates of the text
                    label = "≈ 8y", size = 5)
ggsurv
dev.off()

#creating a table named data with the time in years and with the status that in this case is the occurrence of death
# 1 if it happen and 0 if it didn't happen
time<-apply(ClinicalInfo[,c("days_to_death", "days_to_last_followup")],1 , function(x) as.numeric(na.omit(x)))/365
Status<-as.numeric(ClinicalInfo[, "death_status"])
data<-data.frame("time"=time, "status"=Status)

#creating a survival object and the survival curve
surv_data<-Surv(data$time, data$status)
surv_fit<-survfit(surv_data ~ 1, data = data)
summary(surv_fit)
survfit(formula = surv_data ~ 1, data = data)
#A 50% survival probability isn't reached
min(data$time)
#0.260274, 3 months
max(data$time)
#22.43, 22 years

#ploting the survival curve and saving it in a image
png("OS.png", res=1200, width=7000, height=7000)
ggsurv <- ggsurvplot(surv_fit,
           xlab = "Time (years)",
           ylab = "Survival Probability",
           palette = c("red"),
           conf.int = T,
           pval.coord = c(0,0),
           risk.table = FALSE,
           risk.table.col = " ",
           risk.table.height = 0.25,
           ggtheme = theme_bw(base_family =  "sans"),
           censor = T,
           data=data,
           font.title=c(20, "bold", "black"),
           font.x=c(15, "bold", "black"), 
           font.y=c(15, "bold", "black"), 
           font.tickslab=c(15, "bold", "black"),
           font.legend=c(13, "bold", "black"),
           surv.median.line = "hv"
)

ggsurv
dev.off()


##Survival based on patient's Age groups
#Separating the samples in groups by the median Age
Levels <- ClinicalInfo$Age
min(Levels)
max(Levels)
Group <- ifelse(Levels >= median(Levels), "Higher or equal to 62", "Lower than 62" )
time<-apply(ClinicalInfo[,c("days_to_death", "days_to_last_followup")],1 , function(x) as.numeric(na.omit(x)))/365
Status<-as.numeric(ClinicalInfo[, "death_status"])
data<-data.frame("time"=time, "status"=Status,"group"=Group)

#creating a survival object and the survival curve with the Age groups
surv_data<-Surv(data$time, data$status)
surv_fit_groups<-survfit(surv_data ~ group, data = data)
surv_groups<-survdiff(surv_data ~ group, data=data)
pvalue <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
#p-value=0.0011
summary(surv_fit_groups)
survfit(formula = surv_data ~ group, data = data)
#A 50% survival probability is reached for the ≥ 62y group after 7 years (7y), 
#and not reached for the < 62y

#ploting the survival curve and saving it as an image
png("Age_OS.png", res = 1200, width = 7000, height = 7000)
ggsurv<-ggsurvplot(surv_fit_groups,
           xlab = "Time (years)",
           ylab = "Survival Probability",
           legend.labs=c("≥ 62y", "< 62y"),
           palette = c("red", "blue"),
           conf.int = T,
           conf.int.alpha = 0.05,
           pval = F,
           pval.coord = c(0,0),
           risk.table = FALSE,
           risk.table.col = " ",
           risk.table.height = 0.25,
           ggtheme = theme_bw(base_family =  "sans"),
           censor = T,
           data=data,
           font.title=c(20, "bold", "black"),
           font.x=c(15, "bold", "black"), 
           font.y=c(15, "bold", "black"), 
           font.tickslab=c(15, "bold", "black"),
           font.legend=c(13, "bold", "black"),
           surv.median.line = "hv"
)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 8.5, y = 0.5, # x and y coordinates of the text
                    label = "≈ 7y", size = 5)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate ("text",
                     x = 2, y = 0.01, 
                     label = "p = 0.0011", cex=5, col="black", fontface=3)
ggsurv
dev.off()

##Survival based on tumour Size groups
#Separating the samples in groups by the median Size
Levels <- ClinicalInfo$Size
min(Levels)
max(Levels)
Group <- ifelse(Levels >= median(Levels) , "Greater or equal to 25cm", "Smaller than 25cm")
time<-apply(ClinicalInfo[,c("days_to_death", "days_to_last_followup")],1 , function(x) as.numeric(na.omit(x)))/365
Status<-as.numeric(ClinicalInfo[, "death_status"])
data<-data.frame("time"=time, "status"=Status,"group"=Group)

#creating a survival object and the survival curve with the Size groups
surv_data<-Surv(data$time, data$status)
surv_fit_groups<-survfit(surv_data ~ group, data = data)
surv_groups<-survdiff(surv_data ~ group, data=data)
pvalue <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
#p-value=0.034
summary(surv_fit_groups)
survfit(formula = surv_data ~ group, data = data)
#A 50% survival probability is reached for the ≥ 25cm only after 8.01 (8y) 
#and not reached for the < 25cm 

png("Size_OS.png", res = 1200, width = 7000, height = 7000)
#grDevices::windows()
ggsurv <- ggsurvplot(surv_fit_groups,
           xlab = "Time (years)",
           ylab = "Survival Probability",
           legend.labs=c("≥ 25mm", "< 25mm"),
           palette = c("red", "blue"),
           conf.int = T,
           conf.int.alpha = 0.05,
           pval = F,
           pval.coord = c(0,0),
           risk.table = FALSE,
           risk.table.col = " ",
           risk.table.height = 0.25,
           ggtheme = theme_bw(base_family =  "sans"),
           censor = T,
           data=data,
           font.title=c(20, "bold", "black"),
           font.x=c(15, "bold", "black"), 
           font.y=c(15, "bold", "black"), 
           font.tickslab=c(15, "bold", "black"),
           font.legend=c(13, "bold", "black"),
           font.pval=c(14, "italic", "black"),
           surv.median.line = "hv"
)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 10, y = 0.5, # x and y coordinates of the text
                    label = "≈ 8y", size = 5)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate ("text",
                     x = 1.75, y = 0.01, 
                     label = "p = 0.034", cex=5, col="black", fontface=3)

ggsurv
dev.off()


##Survival based on the presence of metastasis or not groups
#Separating the samples in groups by if there is metastasis or not
MetLevels <- ClinicalInfo$Dist_Met_Dt
MetStatus <- ifelse(is.na(MetLevels), "Without metastasis", "With metastasis")
time<-apply(ClinicalInfo[,c("days_to_death", "days_to_last_followup")],1 , function(x) as.numeric(na.omit(x)))/365
vitalStatus<-as.numeric(ClinicalInfo[, "death_status"])
data<-data.frame("time"=time, "status"=vitalStatus,"group"=MetStatus)

#creating a survival object and the survival curve with the metastasis groups
surv_data<-Surv(data$time, data$status)
surv_fit_groups<-survfit(surv_data ~ group, data = data)
surv_groups<-survdiff(surv_data ~ group, data=data)
pvalue <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
#p<0.0001
summary(surv_fit_groups)
survfit(formula = surv_data ~ group, data = data)
#A 50% survival probability is reached for the group with distant metastasis after 3.9 years (3y), 
#and not reached for the group without metastasis

png("Metastasis_OS.png", res = 1200, width = 7000, height = 7000)
ggsurv<-ggsurvplot(surv_fit_groups,
           xlab = "Time (years)",
           ylab = "Survival Probability",
           legend.labs=c("w/ Metastasis", "w/o Metastasis"),
           palette = c("red", "blue"),
           conf.int = T,
           conf.int.alpha = 0.05,
           pval = F,
           pval.coord = c(0,0),
           risk.table = FALSE,
           risk.table.col = " ",
           risk.table.height = 0.25,
           ggtheme = theme_bw(base_family =  "sans"),
           censor = T,
           data=data,
           font.title=c(20, "bold", "black"),
           font.x=c(15, "bold", "black"), 
           font.y=c(15, "bold", "black"), 
           font.tickslab=c(15, "bold", "black"),
           font.legend=c(13, "bold", "black"),
           surv.median.line = "hv"
)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 5, y = 0.5, # x and y coordinates of the text
                    label = "≈ 4y", size = 5)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate ("text",
                     x = 2, y = 0.01, 
                     label = "p < 0.0001", cex=5, col="black", fontface=3)
ggsurv
dev.off()


##Progression Free Survival based on patient's age group
#Separating the samples in groups by the median Age
AgeLevels <- ClinicalInfo$Age
AgeStatus <- ifelse(AgeLevels >= median(AgeLevels), "Higher or equal to 62", "Lower than 62" )
time <- apply(ClinicalInfo[c("PFS")], 1, function(x) as.numeric(na.omit(x)))/365
vitalStatus <- as.numeric(unlist(ClinicalInfo[,"event"]))
data <- data.frame("time"=time, "status"=vitalStatus, "group"=AgeStatus)
colnames(data)<-c('time','status',"group")

#creating a survival object and the survival curve with the Age groups
surv_data <- Surv(data$time, data$status)
surv_fit_groups <- survfit(surv_data ~ group, data=data)
surv_groups <- survdiff(surv_data ~ group, data=data)
pvalue <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
#p=0.0087
summary(surv_fit_groups)
survfit(formula = surv_data ~ group, data = data)
#A 50% of Progression free survival probability is reached for the group older than 62 after 6.6 years (6y), 
#and not reached for the group younger than 62

png("Age_PFS.png", res = 1200, width = 7000, height = 7000)
ggsurv<-ggsurvplot(surv_fit_groups,
           xlab = "Time (years)",
           ylab = "Progression-Free Survival Probability",
           legend.labs=c("≥ 62y", "< 62y"),
           palette = c("red", "blue"),
           conf.int = T,
           conf.int.alpha = 0.05,
           pval = F,
           pval.coord = c(0,0),
           risk.table = FALSE,
           risk.table.col = " ",
           risk.table.height = 0.25,
           ggtheme = theme_bw(base_family =  "sans"),
           censor = T,
           data=data,
           font.title=c(20, "bold", "black"),
           font.x=c(15, "bold", "black"), 
           font.y=c(15, "bold", "black"), 
           font.tickslab=c(15, "bold", "black"),
           font.legend=c(13, "bold", "black"),
           surv.median.line = "hv"
)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 7.5, y = 0.5, # x and y coordinates of the text
                    label = "≈ 6y", size = 5)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate ("text",
                     x = 1.5, y = 0.01, 
                     label = "p = 0.0087", cex=5, col="black", fontface=3)
ggsurv
dev.off()

##Progression Free Survival based on tumour Size groups
#Separating the samples in groups by the median Size
SizeLevels <- ClinicalInfo$Size
SizeStatus <- ifelse(SizeLevels >= median(SizeLevels), "Greater or equal to 25cm", "Smaller than 25cm")
time <- apply(ClinicalInfo[c("PFS")], 1, function(x) as.numeric(na.omit(x)))/365
vitalStatus <- as.numeric(unlist(ClinicalInfo[,"event"]))
data <- data.frame("time"=time, "status"=vitalStatus, "group"=SizeStatus)
colnames(data)<-c('time','status',"group")

#creating a survival object and the survival curve with the Size groups
surv_data <- Surv(data$time, data$status)
surv_fit_groups <- survfit(surv_data ~ group, data=data)
surv_groups <- survdiff(surv_data ~ group, data=data)
pvalue <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
#p=0.049
summary(surv_fit_groups)
survfit(formula = surv_data ~ group, data = data)
#A 50% of Progression free survival probability is reached for the group smaller than 25cm after 14 years,
#and after 6.68 yearsfor the group bigger or equal to 25cm

png("Size_PFS.png", res = 1200, width = 7000, height = 7000)
ggsurv<-ggsurvplot(surv_fit_groups,
           xlab = "Time (years)",
           ylab = "Progression-Free Survival Probability",
           xlim = c(0,18),
           legend.labs=c("≥ 25mm", "< 25mm"),
           palette = c("red", "blue"),
           conf.int = T,
           conf.int.alpha = 0.05,
           pval = F,
           pval.coord = c(0,0),
           risk.table = FALSE,
           risk.table.col = " ",
           risk.table.height = 0.25,
           ggtheme = theme_bw(base_family =  "sans"),
           censor = T,
           data=data,
           font.title=c(20, "bold", "black"),
           font.x=c(15, "bold", "black"), 
           font.y=c(15, "bold", "black"), 
           font.tickslab=c(15, "bold", "black"),
           font.legend=c(13, "bold", "black"),
           surv.median.line = "hv"
)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 8, y = 0.53, # x and y coordinates of the text
                    label = "≈ 6y", size = 5)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 15.3, y = 0.5, # x and y coordinates of the text
                    label = "≈ 14y", size = 5)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate ("text",
                     x = 1.55, y = 0.01, 
                     label = "p = 0.049", cex=5, col="black", fontface=3)
ggsurv
dev.off()


############################################################## END ##########################################################
