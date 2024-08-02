#install.packages("limma")
#install.packages("openxlsx")
#install.packages("tidyverse")
#install.packages("ggrepel")

#Select the Directory where you can find the RData with the information regarding the TCGA-TNBC cohort
#obtained from the "Select the TNBC.R" code

##############################################################################################
######### Obtain the DEG from the TCGA-TNBC cohort based on ST6GALNAC1 expression#############
##############################################################################################

library(limma)

load("data/Reads_TNBC.RData")

geneExpLevels <- Reads_TNBC["ST6GALNAC1",]
geneExpStatus <- as.data.frame(ifelse(geneExpLevels > median(geneExpLevels), "High", "Low"))
colnames(geneExpStatus) <- "ExpStatus"
Low_ID <- rownames(geneExpStatus)[geneExpStatus$ExpStatus == "Low"]
High_ID <- rownames(geneExpStatus)[geneExpStatus$ExpStatus == "High"]

order <- c(Low_ID, High_ID)
Reads_TNBC_o <- Reads_TNBC[,match(order,colnames(Reads_TNBC))]
refinedReads_TNBC_o <- Reads_TNBC_o[rownames(Reads_TNBC_o) != "" | rownames(Reads_TNBC_o) != NA,]
refinedReads_TNBC_o <- Reads_TNBC_o[rownames(Reads_TNBC_o) != "",]
refinedReads_TNBC_o <- na.omit(refinedReads_TNBC_o)

#Filter low expressed genes
data <- refinedReads_TNBC_o[rowSums(refinedReads_TNBC_o) > 10,]

#Identify the differentially expressed genes between low and high samples
library(openxlsx)
library(tidyverse)
library(ggrepel)

#Create design matrix with coefficients to be used in linear models
sampleType <- as.vector(c(rep("Low", length(Low_ID)), rep("High", length(High_ID))))
design <- model.matrix(~ 0 + sampleType); 
colnames(design) <- gsub("sampleType", "", colnames(design))
rownames(design) <- colnames(data)

#Normalize data
data_norm <- voom(refinedReads_TNBC_o, design)

#Fit a linear model per gene
data_model <- lmFit(data_norm, design)

#Calculate the statistics for our specific comparisons of interest
comparisons <- makeContrasts(High-Low, levels=design)
data_comparisons <- contrasts.fit(data_model, contrasts = comparisons)

#Calculate levels of significance. Uses an empirical Bayes algorithm
#to shrink the gene-specific variances towards the average variance
#across all genes.
data_final<- eBayes(data_comparisons)

#Get differentially expressed genes for each comparison (Moderated T-test)
DEG_TCGA<- topTable(data_final, coef=1,adjust.method="fdr", p.value=1, n="inf")

#Get differentially expressed genes (apply logFC cut-off)
DEG_TCGA <- topTable(data_final, n ="inf", adjust.method="fdr", p.value = 0.05, lfc = 1)
#write.xlsx(DEG_TCGA, "DEG_TCGA_TNBC.xlsx")

allGenes <- rownames(data_final$coefficients)

DEG_TCGA<- topTable(data_final, coef = 1, adjust.method = "fdr", number = "inf")
DEG_TCGA$diffexpressed <- "NO"
DEG_TCGA$diffexpressed[DEG_TCGA$logFC > 1 & DEG_TCGA$P.Value < 0.05] <- "UP"
DEG_TCGA$diffexpressed[DEG_TCGA$logFC < -1 & DEG_TCGA$P.Value < 0.05] <- "DOWN"
head(DEG_TCGA[order(DEG_TCGA$P.Value) & DEG_TCGA$diffexpressed == 'DOWN', ])
top_genes <- subset(DEG_TCGA, (logFC > 1 | logFC < -1) & P.Value < 0.05)

#Sort by absolute logFC descending to select top 20
top_genes <- top_genes[order(-abs(top_genes$logFC)), ]  
top_genes <- top_genes$ID[1:20]

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), 
                                          size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), 
                                          size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
myvolcanoplot <- ggplot(data = DEG_TCGA, aes(x = logFC, y = -log10(P.Value), col = diffexpressed, label = ID)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#2166AC", "black", "#D6604D"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) + 
  coord_cartesian(ylim = c(0, 26), xlim = c(-4, 4)) + 
  labs(color = 'Legend',
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-4, 4, 2)) +
  geom_text_repel(data = subset(DEG_TCGA, ID %in% top_genes),colour = "black" , max.overlaps = Inf) 

png("Volcano Plot.png", res=600, width=5000, height=5000)
myvolcanoplot
dev.off()

############################################################## END ##########################################