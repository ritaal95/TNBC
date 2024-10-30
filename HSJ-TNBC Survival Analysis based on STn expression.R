#############################################################################################################################
############################# Survival Analysis of TNBC Patients Samples Based on STn Expression ############################
#############################################################################################################################

##Libraries used
#install.packages(survival)
#install.packages(survminer)
#install.packages(readxl)

##Load libraries
library(survival)
library(survminer)
library(readxl)

#Select the Directory where you can find the Table S2 with the information regarding the HSJ-TNBC cohort

##Preparing the data table for analysis
TNBC_Data <- as.data.frame(read_xlsx("data/HSJ-TNBC_Data.xlsx"))
ClinicalInfo <-  dplyr::select(TNBC_Data, c("Dx_Date", "Dist_Met_Dt", "Relapse", 
                                            "Relapse_Dt", "Alive", "Death", "Sialyl-Tn_CM_Tumor" ))
rownames(ClinicalInfo) <- TNBC_Data$ID
colnames(ClinicalInfo) <- c("Dx_Date", "Dist_Met_Dt", "Relapse", 
                            "Relapse_Dt", "Alive", "Death", "STn_Score")

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

MarkerExpLevels <- ClinicalInfo$STn_Score
MarkerExpStatus <- ifelse(MarkerExpLevels < "1", "STn-", "STn+")
time<-apply(ClinicalInfo[,c("days_to_death", "days_to_last_followup")],1 , function(x) as.numeric(na.omit(x)))/365
vitalStatus<-as.numeric(ClinicalInfo[, "death_status"])
data<-data.frame("time"=time, "status"=vitalStatus,"group"=MarkerExpStatus)
surv_data<-Surv(data$time, data$status)
surv_fit_groups<-survfit(surv_data ~ group, data = data)
surv_groups<-survdiff(surv_data ~ group, data=data)
pvalue <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
#p-value=0.042
summary(surv_fit_groups)
survfit(formula = surv_data ~ group, data = data)
#A 50% survival probability is reached for the STn+ group after 7.38 years, and not reached for the STn-

png("STn_OS.png", res = 1200, width = 7000, height = 7000)
ggsurv<-ggsurvplot(surv_fit_groups,
           xlab = "Time (years)",
           legend.labs=c("STn-", "STn+"),
           palette = c("blue", "red"),
           conf.int = T,
           conf.int.alpha = 0.05,
           pval = TRUE,
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
                    x = 6.6, y = 0.46, # x and y coordinates of the text
                    label = "≈ 7 y", size = 5)

ggsurv
dev.off()

PFS <- as.data.frame(ifelse(is.na(pmin(ClinicalInfo[,2],ClinicalInfo[,4],ClinicalInfo[,6], na.rm = TRUE)-ClinicalInfo[,1]), ClinicalInfo[,10], pmin(ClinicalInfo[,2],ClinicalInfo[,4], ClinicalInfo[,6],na.rm = TRUE)-ClinicalInfo[,1]))
colnames(PFS) <- "PFS"
event <- as.data.frame(ifelse(is.na(ClinicalInfo[,2]), ifelse(is.na(ClinicalInfo[,4]), ifelse(is.na(ClinicalInfo[,6]), 0, 1), 1), 1))
colnames(event) <- "event"
ClinicalInfo <- cbind(ClinicalInfo, PFS, event)

MarkerExpLevels <- ClinicalInfo$STn_Score
MarkerExpStatus <- ifelse(MarkerExpLevels >= "1", "STn+", "STn-")

time <- apply(ClinicalInfo[c("PFS")], 1, function(x) as.numeric(na.omit(x)))/365
vitalStatus <- as.numeric(unlist(ClinicalInfo[,"event"]))
data <- data.frame("time"=time, "status"=vitalStatus, "group"=MarkerExpStatus)
colnames(data)<-c('time','status',"group")
surv_data <- Surv(data$time, data$status)
surv_fit_groups <- survfit(surv_data ~ group, data=data)
surv_groups <- survdiff(surv_data ~ group, data=data)
pvalue <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
#pValue=0.068
summary(surv_fit_groups)
survfit(formula = surv_data ~ group, data = data)
#50% probability of progression free survival is reached for STn+ group after approximately 5 years, and reached after 14 years for the STn- group

png("STn_PFS.png", res = 1200, width = 7000, height = 7000)
ggsurv<-ggsurvplot(surv_fit_groups,
                   xlab = "Time (years)",
                   ylab = "Progression Free Survival Probability",
                   legend.labs=c("STn-", "STn+"),
                   palette = c("blue", "red"),
                   conf.int = T,
                   conf.int.alpha = 0.05,
                   pval = TRUE,
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
                    x = 4.5, y = 0.46, # x and y coordinates of the text
                    label = "≈ 5 y", size = 5)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 13, y = 0.46, # x and y coordinates of the text
                    label = "≈ 14 y", size = 5)
ggsurv
dev.off()

############################################################## END ##########################################################
