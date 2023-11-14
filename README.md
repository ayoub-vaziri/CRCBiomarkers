# CRCBiomarkers

## Overview
This github repository contains the data files and analysis code used for the scientific paper titled **"Identifying novel gene biomarkers for colorectal cancer diagnosis through bioinformatics analysis and machine learning"**.
The files are organised into four folders:
 - *Data*: which contains all the transcriptomic data required to perform the analyses described in the paper
 - *Codes*: contains the R code to reproduce all  analyses.
 - *Results*: contains all the results produced by the R scripts.
 - *CRCTest*: contains a trained random forest model and code for diagnosing tumor samples from normal samples.

## Reproducing the results
This repository contains all the code necessary to reproduce the results in the paper. In particular:
 - **Codes/DataProcessing** contains the necessary code for preprocessing the initial dataset.
   - Run *mergeDatasets.R* to merge the five microarray datasets (GSE10950, GSE25070, GSE41328, GSE74602, and GSE142279) based on their common genes.
   - Run *batchCorrection.R* to correct batch effects in merged datasets.
   - Run *trainTestSplit.R* to split merged datasets into training and testing datasets.
   
 - **Codes/KeyGenes** includes the required code for identifying key genes.
   - Run *differentialExpressionAnalysis.R* to identify differentially expressed genes (DEGs) on the training dataset.
   - Run *enrichmentAnalysis.R* to gain a deeper insight into the biological significance of the DEGs.
   - Run *geneCoexpressionAnalysis.R* to construct co-expression modules on the training dataset using the automatic network construction package CEMiTool with default settings.
   - Run *overlappedGenes.R* to identify genes overlapped between DEGs and genes within the most significant module identified by CEMiTool.
   - Run *centralityAnalysis.R* to identify key genes in the protein-protein interaction (PPI) network.
   
 - **Codes/DiagnosticGenes** contains the necessary code for identifying and validating diagnostic genes in CRC.
   - Run *LASSO.R* to identify candidate diagnostic genes from within the set of key genes on the training dataset.
   - Run *ROC.R* to evaluate the sensitivity and specificity of candidate diagnostic genes on both the training and testing datasets.
   - Run *boxplot.R* to validate the expression of diagnostic genes on the tesing dataset.
   - Run *trainEvaluateML.R* to verify the accuracy of diagnostic genes in distinguishing between normal and tumor samples using two machine learning models, namely, RF and SVM.
   - Run *externalTesting.R* to ensure the accuracy of diagnostic genes using two external datasets (GSE21815 and GSE106582).
 
## Required software
The scripts use core R functionality and several publicly available R packages listed below. Version numbers in brackets correspond to the versions of the packages that were used to develop and debug these scripts.

 - **R** (4.3.1)
 - **RStudio** (2022.07.2+576): Optional, testing functions and running the code step-by-step.
 - **GEOquery** (2.68.0)
 - **org.Hs.eg.db** (3.17.0)
 - **EnsDb.Hsapiens.v79** (2.99.0)
 - **ggplot2** (3.4.2)
 - **sva** (3.48.0)
 - **igraph** (1.5.0)
 - **CINNA** (1.2.0)
 - **ggraph** (2.1.0)
 - **limma** (3.56.2)
 - **clusterProfiler** (4.8.2)
 - **CEMiTool** (1.24.0)
 - **pheatmap** (1.0.12)
 - **glmnet** (4.1.7)
 - **pROC** (1.18.2)
 - **MLmetrics** (1.1.1)
 - **randomForest** (4.7.1.1)
 - **caret** (6.0.94)
 - **e1071** (1.7.13)
 - **data.table** (1.14.8)
 - **RColorBrewer** (RColorBrewer)
 - **preprocessCore** (preprocessCore)
 - **dplyr** (dplyr)
 - **annotate** (1.78.0)
 - **forcats** (1.0.0)
 - **ggpubr** (0.6.0)
 - **readr** (2.1.4)
 - **graphlayouts** (1.0.0)
 - **enrichplot** (1.20.0)
 - **AnnotationDbi** (1.62.2)
 - **ggvenn** (0.1.10)
 - **ggrepel** (0.9.3)