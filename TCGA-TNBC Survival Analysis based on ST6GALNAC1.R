#install.packages("survival")
#install.packages("survminer")

#Select the Directory where you can find the RData with the information regarding the TCGA-TNBC cohort
#obtained from the "Select the TNBC.R" code

#############################################################################################################
########################### Complete Overall Survival Analysis of TNBC TCGA Cohort ##########################
#############################################################################################################

load("data/TNBC_data.RData")

library(survival)
library(survminer)


tnbc.vars <- c(
  "year_of_initial_pathologic_diagnosis",
  'vital_status','days_to_death','days_to_last_followup',
  'age_at_initial_pathologic_diagnosis','tissue_source_site',
  "tumor_tissue_site",'primary_lymph_node_presentation_assessment',
  'lymph_node_examined_count',"anatomic_neoplasm_subdivisions",
  'number_of_lymphnodes_positive_by_ihc','number_of_lymphnodes_positive_by_he',
  'axillary_lymph_node_stage_method_type',
  'axillary_lymph_node_stage_other_method_descriptive_text',
  'histological_type','histological_type_other',
  'cytokeratin_immunohistochemistry_staining_method_micrometastasis_indicator',
  "margin_status",'first_nonlymph_node_metastasis_anatomic_sites',
  'distant_metastasis_present_ind2','radiation_therapy',
  'stage_event_pathologic_stage','stage_event_tnm_categories',
  "stage_event_clinical_stage",'surv_event_time',"gender"
)


clinicalInfo <- clinical$primary.solid.tumor[, tnbc.vars]

patientsIDs <- colnames(logTPMs)

time <- apply(clinicalInfo[patientsIDs, c("days_to_death","days_to_last_followup")], 1, function(x) as.numeric(na.omit(x)))/365
vitalStatus <- as.numeric(clinicalInfo[patientsIDs,"vital_status"])
data<-data.frame("time"=time, "status"=vitalStatus)
surv_data<-Surv(data$time, data$status)
surv_fit<-survfit(surv_data ~ 1, data = data)
summary(surv_fit)
survfit(formula = surv_data ~ 1, data = data)
#A 50% survival probability is reached after 8 years
median(data$time)
#1.293151, 1 year and approximately 4 months of median follow-up time
min(data$time)
#0
max(data$time)
# 9.512329

png("TCGA_TNBC_OS.png", res=1200, width=7000, height=7000)
#grDevices::windows()
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
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 8.7, y = 0.5, # x and y coordinates of the text
                    label = "≈ 8 y", size = 5)
ggsurv
dev.off()

#Survival by ST6GALNAC1 high and low
geneExpLevels <- logTPMs["ST6GALNAC1",]
geneExpStatus <- ifelse(geneExpLevels > median(geneExpLevels), "High", "Low")

clinicalInfo$vital_status <- ifelse(clinicalInfo$vital_status == "Alive", 0, 1)

time <- apply(clinicalInfo[patientsIDs, c("days_to_death","days_to_last_followup")], 1, function(x) as.numeric(na.omit(x)))/365
vitalStatus <- as.numeric(clinicalInfo[patientsIDs,"vital_status"])
data <- data.frame("time"=time, "status"=vitalStatus, "group"=geneExpStatus)
surv_data <- Surv(data$time, data$status)
surv_fit_groups <- survfit(surv_data ~ group, data=data)
surv_groups <- survdiff(surv_data ~ group, data=data)
pvalue <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
#pvalue=0.9051407
summary(surv_fit_groups)
survfit(formula = surv_data ~ group, data = data)
#A 50% of disease free survival probability is reached for the low expression group after 8.1 years.


png("ST6GALNAC1_OS.png", res = 1200, width = 7000, height = 7000)
#grDevices::windows()
ggsurv<-ggsurvplot(surv_fit_groups,
                   xlab = "Time (years)",
                   ylab = "Survival Probability",
                   legend.labs = c("ST6GALNAC1 High", "ST6GALNAC1 Low"),
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
                    x = 8.7, y = 0.47, # x and y coordinates of the text
                    label = "≈ 8 y", size = 5)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate ("text",
                     x = 0.8, y = 0.01, 
                     label = "p = 0.91", cex=5, col="black", fontface=3)
ggsurv
dev.off()

############################################################## END ##########################################