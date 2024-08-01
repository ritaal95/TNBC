#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("biomaRt")
#install.packages("openslsx")

#Select the Directory where you can find the RData with the information retrieved from the TCGA
#obtained in the "Dowload TCGA-BRCA.R" code

#############################################################################################################
########################### Selecting the TNBC cases from the TCGA-BRCA project #############################
#############################################################################################################


#Load libraries
library(biomaRt)
library(openxlsx)

#Load data 

load("data/clinical.rda")
load("data/tpm.rda")
load("data/fpkm.per.tissue.rda")
load("data/reads.rda")

#Construction of TNBC dataset
tnbc.vars <- c('breast_carcinoma_estrogen_receptor_status',
               'breast_carcinoma_progesterone_receptor_status',
               'her2_immunohistochemistry_level_result',
               'lab_proc_her2_neu_immunohistochemistry_receptor_status',
               'lab_procedure_her2_neu_in_situ_hybrid_outcome_type'
)

clinical.tnbc <- clinical$primary.solid.tumor[,tnbc.vars]

clinical.tnbc.aux <- clinical.tnbc
names(clinical.tnbc.aux) <- c("ESR1","PGR","HER2_level","HER2_status","HER2_FISH")

summary(clinical.tnbc.aux)
clinical.tnbc.aux$ESR1[clinical.tnbc.aux$ESR1==""] <- NA
clinical.tnbc.aux$PGR[clinical.tnbc.aux$PGR==""] <- NA
clinical.tnbc.aux$HER2_level[clinical.tnbc.aux$HER2_level==""] <- NA
clinical.tnbc.aux$HER2_status[clinical.tnbc.aux$HER2_status==""] <- NA
clinical.tnbc.aux$HER2_FISH[clinical.tnbc.aux$HER2_FISH==""] <- NA

#delete samples with NA expression of ESR
clinical.tnbc.aux <- clinical.tnbc.aux[-(which(is.na(clinical.tnbc.aux[,1]))),] 

#Finding the concordance between the HER2 variables measured.
#Concordance between HER2 (IHC) level and HER2 (IHC) status:

table(as.factor(clinical.tnbc.aux[,3]),clinical.tnbc.aux[,4])

#Creating a vector of HER2 (IHC) status based on HER2 (ICH) level:

clinical.tnbc.aux$my_HER2_level <- as.character(clinical.tnbc.aux$HER2_level)

clinical.tnbc.aux[which(clinical.tnbc.aux$HER2_level == '0'),6] <- 'Negative'
clinical.tnbc.aux[which(clinical.tnbc.aux$HER2_level == '1+'),6] <- 'Negative'
clinical.tnbc.aux[which(clinical.tnbc.aux$HER2_level == '2+'),6] <- 'Indeterminate'
clinical.tnbc.aux[which(clinical.tnbc.aux$HER2_level == '3+'),6] <- 'Positive'

#Individuals for which HER2 (IHC) status different from myHER2 (IHC) status (based on HER2 IHC level)

HER2status.dif.HER2level <- clinical.tnbc.aux[c(rownames(clinical.tnbc.aux[which(clinical.tnbc.aux[,6] == 'Negative' & clinical.tnbc.aux[,4] == 'Positive'),]),
                                                rownames(clinical.tnbc.aux[which(clinical.tnbc.aux[,6] == 'Positive' & clinical.tnbc.aux[,4] == 'Negative'),])),]

# getting the individuals short name, as in the clinical data
clinical.tnbc.barcode.pattern <- '(TCGA-[A-Z0-9a-z]{2}-[a-zA-Z0-9]{4})-([0-9]{2}).+'
clinical.tnbc.shortname <- gsub(clinical.tnbc.barcode.pattern,'\\1',rownames(clinical.tnbc))

# individuals with discordant HER2 (IHC) status and HER2 (IHC) level
HER2status.dif.HER2level.names <- clinical.tnbc.shortname[which(clinical.tnbc.shortname %in% rownames(HER2status.dif.HER2level))]
length(HER2status.dif.HER2level.names)

#Concordance between HER2 (IHC) status and HER2 (FISH):
  
table(clinical.tnbc.aux[,5],clinical.tnbc.aux[,4])

#Individuals for which HER2 (IHC) status different from HER2 (FISH):

HER2status.dif.FISH <- clinical.tnbc.aux[c(rownames(clinical.tnbc.aux[which(clinical.tnbc.aux[,4] == 'Negative' & clinical.tnbc.aux[,5] == 'Positive'),]),
                                           rownames(clinical.tnbc.aux[which(clinical.tnbc.aux[,4] == 'Positive' & clinical.tnbc.aux[,5] == 'Negative'),])),]

HER2status.dif.HER2FISH.names <- clinical.tnbc.shortname[which(clinical.tnbc.shortname %in% rownames( HER2status.dif.FISH))]
length(HER2status.dif.HER2FISH.names)

#suspect individuals with discordante HER2 labels from different lab testing
suspect <- c(HER2status.dif.HER2level.names,HER2status.dif.HER2FISH.names)
length(suspect)

#Replacing HER2 (IHC) status by HER2 (FISH) (a more accurate test) whenever available:
  
clinical.tnbc.aux$HER2_status_plus_FISH <- clinical.tnbc.aux[,4]
clinical.tnbc.aux[!is.na(clinical.tnbc.aux[,5]) ,7] <- clinical.tnbc.aux[!is.na(clinical.tnbc.aux[,5]),5]

#Concordance between HER2 (IHC) status and HER2 (IHC status + FISH):
  
table(clinical.tnbc.aux[,4],clinical.tnbc.aux[,7])

#Building the final clinical TNBC dataset, composed of all individuals and gene expression status
#of ESR1, PGR and HER2:
  
clinical.tnbc <- clinical.tnbc.aux[,c(1,2,7)]

# replace all '' or NA with 'Indeterminate'
clinical.tnbc[is.na(clinical.tnbc) | clinical.tnbc == '' | clinical.tnbc == 'Equivocal'] <- 'Indeterminate'

# individuals with at least one 'Positive' and keep all are marked as not TNBC
tnbc.status.pos <- clinical.tnbc == 'Positive'

#tnbc.status.pos <- clinical.tnbc[,c(1,2,3)] == 'Positive'
not.tnbc.id <- rownames(tnbc.status.pos[which(sapply(seq(nrow(tnbc.status.pos)), function(ix) { any(tnbc.status.pos[ix,]) })),])

# individuals with at least one 'Indeterminate' and excluded from dataset
tnbc.status.ind <- clinical.tnbc[!(rownames(clinical.tnbc) %in% not.tnbc.id),] == 'Indeterminate'
not.tnbc.id     <- rownames(tnbc.status.ind[which(sapply(seq(nrow(tnbc.status.ind)), function(ix) { any(tnbc.status.ind[ix,]) })),])
clinical.tnbc <- clinical.tnbc[!(rownames(clinical.tnbc) %in% not.tnbc.id), ]

tnbc.status.neg <- clinical.tnbc == 'Negative'
tnbc.ix <- sapply(seq(nrow(clinical.tnbc)), function(ix) { all(tnbc.status.neg[ix,]) })

# set cases with TNBC
clinical.tnbc$tnbc <- 'NO_TNBC'
clinical.tnbc[tnbc.ix, 'tnbc'] <- 'TNBC'

#Building the Y (binary data, i.e., TNBC vs. nonTNBC) response vector:
# building the ydata (binary)
ydata <- data.frame(tnbc = clinical.tnbc$tnbc, row.names = rownames(clinical.tnbc))
ydata$tnbc <- factor(ydata$tnbc)

TNBCsamples<-clinical.tnbc[which(clinical.tnbc$tnbc == "TNBC"),] #160 TNBC samples 
TNBCsamples<-rownames(TNBCsamples)
length(TNBCsamples) 
#160 TNBC samples 

NoTNBCsamples <- clinical.tnbc[which(clinical.tnbc$tnbc == "NO_TNBC"),]
NoTNBCsamples<-rownames(NoTNBCsamples)
length(NoTNBCsamples)
#863 BC samples not TNBC

ClinicalDataTNBC<-clinical.tnbc[TNBCsamples,]

TPM_PrimaryTumor <- tpm$primary.solid.tumor
TPM_SolidNormal <- tpm$solid.tissue.normal

reads_PrimaryTumor <- reads$primary.solid.tumor

ClinicalDataNotTNBC <- clinical.tnbc[NoTNBCsamples,]

###############################
#### Gene Names Annotation ####
###############################

test <- rownames(TPM_PrimaryTumor)
rows <- c()
for(i in 1:length(test)){
  rows <- rbind(rows, unlist(strsplit(as.character(rownames(TPM_PrimaryTumor)[i]), split=".", fixed=T))[1])
}

rownames(TPM_PrimaryTumor)<-rows

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

attributes = listAttributes(ensembl)
#attributes[1:50,]

geneIDs <- rownames(TPM_PrimaryTumor)
length(rownames(TPM_PrimaryTumor))
allInfo <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                 filters = "ensembl_gene_id", 
                 values = geneIDs,
                 mart = ensembl)

geneNames <- allInfo[match(geneIDs, allInfo[,2]),1]
rownames(TPM_PrimaryTumor)<-geneNames

SamplesName <- colnames(TPM_PrimaryTumor)
SampleNameShort<- gsub("-01.*", "", SamplesName)
colnames(TPM_PrimaryTumor)<-SampleNameShort
TPM_PrimaryTumor<-TPM_PrimaryTumor+1
logTPMs_All <- log2(TPM_PrimaryTumor)
logTPMs <- logTPMs_All[, TNBCsamples]

test <- rownames(TPM_SolidNormal)
rows <- c()
for(i in 1:length(test)){
  rows <- rbind(rows, unlist(strsplit(as.character(rownames(TPM_SolidNormal)[i]), split=".", fixed=T))[1])
}

rownames(TPM_SolidNormal)<-rows

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

attributes = listAttributes(ensembl)
#attributes[1:50,]

geneIDs <- rownames(TPM_SolidNormal)
length(rownames(TPM_SolidNormal))
allInfo <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                 filters = "ensembl_gene_id", 
                 values = geneIDs,
                 mart = ensembl)

geneNames <- allInfo[match(geneIDs, allInfo[,2]),1]
rownames(TPM_SolidNormal)<-geneNames

SamplesName <- colnames(TPM_SolidNormal)
SampleNameShort <- gsub("-11.*", "", SamplesName)
colnames(TPM_SolidNormal)<-SampleNameShort
TPM_SolidNormal <- TPM_SolidNormal+1
logTPM_Normal <- log2(TPM_SolidNormal)
logTPM_NormalTNBC <- logTPM_Normal[,match(TNBCsamples, colnames(TPM_SolidNormal))]
logTPM_NormalTNBC <- logTPM_NormalTNBC[ ,!apply(logTPM_NormalTNBC, 2, function(x) all(is.na(x)) ) ]

##TNBC Samples that have normal cases associated
TNBCsamples_Normal <- colnames(logTPM_NormalTNBC)
logTPMs_TNBC_Normal <- logTPMs[,TNBCsamples_Normal]


##BC without TNBC
logTPMs_BC <- logTPMs_All[, NoTNBCsamples]
logTPM_NormalBC <- logTPM_Normal[,match(NoTNBCsamples, colnames(TPM_SolidNormal))]
logTPM_NormalBC <- logTPM_NormalBC[ ,!apply(logTPM_NormalBC, 2, function(x) all(is.na(x)) ) ]

##BC Samples that have normal cases associated
BCsamples_Normal <- colnames(logTPM_NormalBC)
logTPMs_BC_Normal <- logTPMs_BC[,BCsamples_Normal]

##READS
test <- rownames(reads_PrimaryTumor)
rows <- c()
for(i in 1:length(test)){
  rows <- rbind(rows, unlist(strsplit(as.character(rownames(reads_PrimaryTumor)[i]), split=".", fixed=T))[1])
}

rownames(reads_PrimaryTumor)<-rows

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

attributes = listAttributes(ensembl)
#attributes[1:50,]

geneIDs <- rownames(reads_PrimaryTumor)
length(rownames(reads_PrimaryTumor))
allInfo <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id"),
                 filters = "ensembl_gene_id", 
                 values = geneIDs,
                 mart = ensembl)

geneNames <- allInfo[match(geneIDs, allInfo[,2]),1]
rownames(reads_PrimaryTumor)<-geneNames

SamplesName <- colnames(reads_PrimaryTumor)
SampleNameShort<- gsub("-01.*", "", SamplesName)
colnames(reads_PrimaryTumor)<-SampleNameShort
Reads_TNBC <- reads_PrimaryTumor[, TNBCsamples]

ls()
write.xlsx(logTPMs, file = "logTPMs.xlsx")
write.xlsx(Reads_TNBC, file = "Reads_TNBC-TCGA.xlsx")

save(logTPMs_All, file="logTPMs.RData")
save(logTPMs, clinical, file = "TNBC_Data.RData")
save(Reads_TNBC, file = "Reads_TNBC.RData")
