# Sialyl-Tn in triple negative breast cancer: gaining insights into a novel subgroup.
**Autors:** Rita Adubeiro Lourenço, Daniela Barreira, Carla Lopes, Pedro Granjo, Ana Sofia Rodrigues, Zélia Silva, Manuela Martins, Ana Rita Grosso, Paula A Videira

# Table of Contents

- [Project Aim](#project-aim)
- [Requirements](#requirements)
- [Performed Analysis](#performed-analysis-for-the-paper)

# Project Aim

The study aimed to elucidate the expression, clinical implications, and biological context of the sialyl-Tn (STn) antigen in Triple Negative Breast Cancer (TNBC). With that in mind, two types of cohorts and datasets were used in this analysis. The first one, named HSJ-TNBC, comprises a dataset of tumour tissue microarray staining data and clinical information and the second one, named TCGA-TNBC, is a dataset retrieved from The Cancer Genome Atlas (TCGA) database and includes transcriptomic and clinical data.

# Requirements

To run the scripts it is required R (programming language) with the following packages installed:

- survival
- survminer
- readxl
- tidyverse
- openxlsx
- corrplot
- BiocManager
- TCGAbiolinks
- SummarizedExperiment
- futile.logger
- EnsDb.Hsapiens.v86
- biomaRt
- limma
- ggrepel


# Performed Analysis for the paper

- Overall and progression-free survival analysis
- Contingency tables analysis
- Correlation analysis
- Comparison analysis
- Retrieving data from the TCGA database
- Selecting a subgroup of tumours (TNBC)
- Differential Expression Analysis
