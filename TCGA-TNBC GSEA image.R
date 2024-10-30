#############################################################
############## Plot gene ontology enrichment ################
#############################################################

##Libraries used
#install.packages("tidyverse")

#Load librarie
#library(tidyverse)

##Select the Directory where you can find the excel with the information regarding the GO Biological processes obtained

# Import the table containing all the enriched GO terms
GO_all <- read.csv("data/GO_BP.csv",header=T,stringsAsFactors = T)

#Select the 30 most enriched (exclude negative NES)
GO_30 <- GO_all[1:31,]
GO_30 <- GO_30[-1,]
GO_30$geneSet <- with(GO_30, paste0("(", geneSet, ")"))
GO_30 <- unite(GO_30, GO_biological_process,
               description,
               geneSet,
               sep = " ")
GO_30 <- GO_30[-c(2,3,5,7,8,10,11)]
GO_30 <- select(GO_30, GO_biological_process, leadingEdgeNum, normalizedEnrichmentScore, FDR)
colnames(GO_30) <- c("GO_biological_process","Gene_number","NES","FDR")
write.table(GO_30, "data/GO_30.txt", sep = "\t", row.names = FALSE)

# Transform the column 'Gene_number' into a numeric variable
GO_30$Gene_number <- as.numeric(GO_30$Gene_number)

# Transform FDR values by -log10('FDR values')
GO_30$'|log10(FDR)|' <- -(log10(GO_30$FDR))
GO_30$GO_biological_process <- gsub("cellular response to vascular endothelial growth factor stimulus", "cellular response to VEGF stimulus", GO_30$GO_biological_process )
GO_30$GO_biological_process <- gsub("phospholipase C-activating G protein-coupled receptor signaling pathway", "PLC-activating GPCR signaling pathway", GO_30$GO_biological_process )
GO_30$GO_biological_process <- gsub("protein kinase A signaling", "PKA signaling", GO_30$GO_biological_process )

 
# Draw the plot with ggplot2
#--------------------------------------
png("data/GO_BP_30.png", res=800, width=10000, height=6300)
#windows()
ggplot(GO_30, aes(x = GO_biological_process, y = NES)) +
  #geom_hline(yintercept = 1, linetype="dashed", 
  #          color = "azure4", size=.5)+
  geom_point(data=GO_30,aes(x=GO_biological_process, y=NES,size = Gene_number, colour = `|log10(FDR)|`), alpha=.7)+
  geom_density(alpha = 0.2) +
  scale_x_discrete(limits= GO_30$GO_biological_process)+
  scale_color_gradient(low="#2166AC",high="#D6604D",limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt"), color = "black", size = 20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"), color = "black", size = 20),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title.align=0.5)+
  xlab("GO biological processes")+
  ylab("Normalized Enrichment Score")+
  labs(color="-log10(FDR)", size="Number\nof genes")+ #Replace by your variable names; \n allow a new line for text
  guides(y = guide_axis(order=2),
         colour = guide_colourbar(order=1))
dev.off()
