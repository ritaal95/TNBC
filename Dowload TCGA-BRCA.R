#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("TCGAbiolinks")
#BiocManager::install("SummarizedExperiment")
#install.packages("futile.logger")


#############################################################################################################
########## Downloading the gene expression and clinical metadata of the TCGA-BRCA ###########################
#############################################################################################################

library("TCGAbiolinks")
library("SummarizedExperiment")
library(futile.logger)


project <- 'TCGA-BRCA'
cached.file <- 'tcga-biolinks-cache.RData'
if (file.exists(cached.file)) {
  load(cached.file)
} else {
  query <- list()
  
  # Query the GDC to obtain clinical data in BCR XML format
  query$clinical <- GDCquery(project = project,
                             data.category = "Clinical",
                             data.format = "BCR XML" )
  # Download clinical data
  GDCdownload(query$clinical, method = 'api')

  # Prepare clinical data
  gdc <- list()
  gdc$clinical        <- GDCprepare_clinic(query$clinical,clinical.info =  'patient')
  gdc$drug            <- GDCprepare_clinic(query$clinical,clinical.info =  'drug')
  gdc$radiation       <- GDCprepare_clinic(query$clinical,clinical.info =  'radiation')
  gdc$admin           <- GDCprepare_clinic(query$clinical,clinical.info =  'admin')
  gdc$follow.up       <- GDCprepare_clinic(query$clinical,clinical.info =  'follow_up')
  gdc$stage_event     <- GDCprepare_clinic(query$clinical,clinical.info =  'stage_event')
  gdc$new_tumor_event <- GDCprepare_clinic(query$clinical,clinical.info =  'new_tumor_event')
  
  # Query the GDC to obtain RNA-Seq data
  query$rnaseq <- GDCquery(project = project,
                           data.category = "Transcriptome Profiling",
                           data.type     = "Gene Expression Quantification",
                           workflow.type = "STAR - Counts")
  
  #Download RNA-Seq data
  GDCdownload(query$rnaseq, method = 'api')
  
  #Prepare RNA-Seq data
  gdc$rnaseq <- GDCprepare(query$rnaseq)
  
  #Save the data to a file
  save(gdc, query,
       file = cached.file)
}  

#Function to get participant code from full barcode
getParticipantCode <- function(fullBarcode) {
  tcga.barcode.pattern <- '(TCGA-[A-Z0-9a-z]{2}-[a-zA-Z0-9]{4}).*'
  return(gsub(tcga.barcode.pattern, '\\1', fullBarcode))
}

getParticipant <- function(fullBarcode) {
getParticipantCode(fullBarcode)
}

#Function to get sample type code from full barcode
getSampleTypeCode <- function(fullBarcode) {
  tcga.barcode.pattern <- '(TCGA-[A-Z0-9a-z]{2}-[a-zA-Z0-9]{4})-([0-9]{2}).*'
  return(gsub(tcga.barcode.pattern, '\\2', fullBarcode))
}

getSample <- function(fullBarcode) {
  getSampleTypeCode(fullBarcode)
}

#Organize RNA-Seq data by participant barcode and sample type
fpkm.per.tissue.barcode <- list()
fpkm.per.tissue.barcode$all <- getParticipant(gdc$rnaseq@colData$barcode)

fpkm.per.tissue <- list()
fpkm.per.tissue$all <- as.matrix(gdc$rnaseq@assays@data$fpkm_unstrand)
colnames(fpkm.per.tissue$all) <- gdc$rnaseq@colData$barcode
rownames(fpkm.per.tissue$all) <- gdc$rnaseq@rowRanges$gene_id

tpm <- list()
tpm$all <- as.matrix(gdc$rnaseq@assays@data$tpm_unstrand)
colnames(tpm$all) <- gdc$rnaseq@colData$barcode
rownames(tpm$all) <- gdc$rnaseq@rowRanges$gene_id

reads <- list()
reads$all <- as.matrix(gdc$rnaseq@assays@data$unstranded)
colnames(reads$all) <- gdc$rnaseq@colData$barcode
rownames(reads$all) <- gdc$rnaseq@rowRanges$gene_id

gene.ranges <- gdc$rnaseq@rowRanges
gene.names <- gdc$rnaseq@rowRanges$gene_name

#Organize clinical data
clinical <- list()
clinical$all <- gdc$clinical
clinical$all <- clinical$all[!duplicated(clinical$all[,1]), ]
rownames(clinical$all) <- clinical$all$bcr_patient_barcode
clinical$all <- clinical$all[,c(-1)]

FollowUp <- list()
FollowUp$all <- gdc$follow.up
rownames(FollowUp$all) <- FollowUp$all$bcr_followup_barcode
FollowUp$all <- FollowUp$all[,c(-2)]

NewTumorEvent <- list()
NewTumorEvent$all <- gdc$new_tumor_event
NewTumorEvent$all <- NewTumorEvent$all[!duplicated(NewTumorEvent$all[,1]), ]
rownames(NewTumorEvent$all) <- NewTumorEvent$all$bcr_patient_barcode
NewTumorEvent$all <- NewTumorEvent$all[,c(-1)]

Drug<- list()
Drug$all <- gdc$drug
Drug$all <- Drug$all[!duplicated(Drug$all[,1]), ]
rownames(Drug$all) <- Drug$all$bcr_patient_barcode
Drug$all <- Drug$all[,c(-1)]

#Map tissue type codes to names
tissue.type  <- as.numeric(getSample(colnames(fpkm.per.tissue$all)))

tissue.mapping <- list( 
  'primary.solid.tumor' = 01,
  'recurrent.solid.tumor' = 02,
  'primary.blood.derived.cancer.peripheral.blood' = 03,
  'recurrent.blood.derived.cancer.bone.marrow' = 04,
  'additional.new.primary' = 05,
  'metastatic' = 06,
  'additional.metastatic' = 07,
  'human.tumor.original.cells' = 08,
  'primary.blood.derived.cancer.bone.marrow' = 09,
  'blood.derived.normal' = 10,
  'solid.tissue.normal' = 11,
  'buccal.cell.normal' = 12,
  'ebv.immortalized.normal' = 13,
  'bone.marrow.normal' = 14,
  'sample.type.15' = 15,
  'sample.type.16' = 16,
  'control.analyte' = 20,
  'recurrent.blood.derived.cancer.peripheral.blood' = 40,
  'cell.lines' = 50,
  'primary.xenograft.tissue' = 60,
  'cell.line.derived.xenograft.tissue' = 61,
  'sample.type.99' = 99
)

#Assign FPKM data to respective tissue types
tissue.ix <- list()
for (el in names(tissue.mapping)) {
  tissue.ix[[el]] <- (tissue.type == tissue.mapping[[el]])
  if (any(tissue.ix[[el]]))
    fpkm.per.tissue[[el]] <- fpkm.per.tissue$all[,tissue.ix[[el]]]
}

#Calculate sample size for FPKM data
sample.size <- c()
for (el in names(fpkm.per.tissue)) {
  sample.size <- rbind(sample.size, c(ncol(fpkm.per.tissue[[el]]), nrow(fpkm.per.tissue[[el]])))
}
rownames(sample.size) <- names(fpkm.per.tissue)
colnames(sample.size) <- c('# of Samples', '# of Genes')
futile.logger::flog.info('Tissue information per tissue type:', sample.size, capture = TRUE)
#primary.solid.tumor 1106 samples and 60660 genes
#metastatic 7 samples and 60660 genes
#solid.tissue.normal 113 samples and 60660 genes

#Store barcodes per tissue type
for (el in names(fpkm.per.tissue)) {
  fpkm.per.tissue.barcode[[el]] <- getParticipant(colnames(fpkm.per.tissue[[el]]))
}

#Process survival event data
survival.event.tmp <- clinical$all[, c('days_to_death', "days_to_last_followup", 
                                       "days_to_last_known_alive", "vital_status")]
vital.status.ix <- !is.na(clinical$all$days_to_death)
clinical$all$surv_event_time[ vital.status.ix] <- clinical$all[ vital.status.ix,'days_to_death']
clinical$all$surv_event_time[!vital.status.ix] <- clinical$all[!vital.status.ix,'days_to_last_followup']

#Organize clinical data by tissue type
for (el in names(fpkm.per.tissue)[names(fpkm.per.tissue) != 'all']) {
  clinical[[el]] <- clinical$all[unique(fpkm.per.tissue.barcode[[el]]),]
}

sample.size <- c()
for (el in names(clinical)) {
  sample.size <- rbind(sample.size, c(nrow(clinical[[el]]), ncol(clinical[[el]])))
}
rownames(sample.size) <- names(fpkm.per.tissue)
colnames(sample.size) <- c('# of Samples', '# of Features')
futile.logger::flog.info('Clinical information per tissue type:', sample.size, capture = TRUE)
#primary.solid.tumor 1095 samples and 114 features
#metastatic 7 samples and 114 features
#solid.tissue.normal 113 samples and 114 features

tissue.type  <- as.numeric(getSample(colnames(tpm$all)))

tissue.mapping <- list( 
  'primary.solid.tumor' = 01,
  'recurrent.solid.tumor' = 02,
  'primary.blood.derived.cancer.peripheral.blood' = 03,
  'recurrent.blood.derived.cancer.bone.marrow' = 04,
  'additional.new.primary' = 05,
  'metastatic' = 06,
  'additional.metastatic' = 07,
  'human.tumor.original.cells' = 08,
  'primary.blood.derived.cancer.bone.marrow' = 09,
  'blood.derived.normal' = 10,
  'solid.tissue.normal' = 11,
  'buccal.cell.normal' = 12,
  'ebv.immortalized.normal' = 13,
  'bone.marrow.normal' = 14,
  'sample.type.15' = 15,
  'sample.type.16' = 16,
  'control.analyte' = 20,
  'recurrent.blood.derived.cancer.peripheral.blood' = 40,
  'cell.lines' = 50,
  'primary.xenograft.tissue' = 60,
  'cell.line.derived.xenograft.tissue' = 61,
  'sample.type.99' = 99
)
tissue.ix <- list()
for (el in names(tissue.mapping)) {
  tissue.ix[[el]] <- (tissue.type == tissue.mapping[[el]])
  if (any(tissue.ix[[el]]))
    tpm[[el]] <- tpm$all[,tissue.ix[[el]]]
}
sample.size <- c()
for (el in names(tpm)) {
  sample.size <- rbind(sample.size, c(ncol(tpm[[el]]), nrow(tpm[[el]])))
}
rownames(sample.size) <- names(tpm)
colnames(sample.size) <- c('# of Samples', '# of Genes')
futile.logger::flog.info('Tissue information per tissue type:', sample.size, capture = TRUE)
#primary.solid.tumor 1106 samples and 60660 genes
#metastatic 7 samples and 60660 genes
#solid.tissue.normal 113 samples and 60660 genes

tissue.type  <- as.numeric(getSample(colnames(reads$all)))

tissue.mapping <- list( 
  'primary.solid.tumor' = 01,
  'recurrent.solid.tumor' = 02,
  'primary.blood.derived.cancer.peripheral.blood' = 03,
  'recurrent.blood.derived.cancer.bone.marrow' = 04,
  'additional.new.primary' = 05,
  'metastatic' = 06,
  'additional.metastatic' = 07,
  'human.tumor.original.cells' = 08,
  'primary.blood.derived.cancer.bone.marrow' = 09,
  'blood.derived.normal' = 10,
  'solid.tissue.normal' = 11,
  'buccal.cell.normal' = 12,
  'ebv.immortalized.normal' = 13,
  'bone.marrow.normal' = 14,
  'sample.type.15' = 15,
  'sample.type.16' = 16,
  'control.analyte' = 20,
  'recurrent.blood.derived.cancer.peripheral.blood' = 40,
  'cell.lines' = 50,
  'primary.xenograft.tissue' = 60,
  'cell.line.derived.xenograft.tissue' = 61,
  'sample.type.99' = 99
)
tissue.ix <- list()
for (el in names(tissue.mapping)) {
  tissue.ix[[el]] <- (tissue.type == tissue.mapping[[el]])
  if (any(tissue.ix[[el]]))
    reads[[el]] <- reads$all[,tissue.ix[[el]]]
}
sample.size <- c()
for (el in names(reads)) {
  sample.size <- rbind(sample.size, c(ncol(reads[[el]]), nrow(reads[[el]])))
}
rownames(sample.size) <- names(reads)
colnames(sample.size) <- c('# of Samples', '# of Genes')
futile.logger::flog.info('Tissue information per tissue type:', sample.size, capture = TRUE)
#primary.solid.tumor 1106 samples and 60660 genes
#metastatic 7 samples and 60660 genes
#solid.tissue.normal 113 samples and 60660 genes

#save(gdc, file = "gdc.rda")
save(fpkm.per.tissue, file="fpkm.per.tissue.rda")
#save(fpkm.per.tissue.barcode, file="fpkm.per.tissue.barcode.rda")
save(clinical, file="clinical.rda")
#save(gene.ranges, file="gene.ranges.rda")
#save(gene.names, file="gene.names.rda")
save(tpm,file="tpm.rda")
#save(cached.file, file="cached.file.rda")
#save(FollowUp, file="FollowUp.rda")
#save(NewTumorEvent, file="NewTumorEvent.rda")
save(reads, file = "reads.rda")

