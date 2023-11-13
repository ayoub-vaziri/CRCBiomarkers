# CRCBiomarkers

## Overview
This github repository contains the data files and analysis code used for the scientific paper titled **"Identifying novel gene biomarkers for colorectal cancer diagnosis through bioinformatics analysis and machine learning"**.
The files are organised into four folders:
 - *Data*: which contains all the transcriptomic data required to perform the analyses described in the paper
 - *Codes*: contains the R code to reproduce all  analyses.
 - *Results*: contains all the results produced by the R scripts.
 - *CRCTest*: contains a trained random forest model and code for diagnosing tumor samples from normal samples.

## Reproducing the results
This repository contains all the code necessary to reproduce the results in the paper:
 - **Codes/DataProcessing** contains the necessary code for preprocessing the initial dataset.
   - Run *mergeDatasets.R* to merge the five microarray datasets (GSE10950, GSE25070, GSE41328, GSE74602, and GSE142279) based on their common genes.
   - Run *batchCorrection.R* to correct batch effects in merged datasets.
   - Run *trainTestSplit.R* to split merged datasets into training and testing datasets.
   
 - **Codes/KeyGenes** includes the required code for identifying key genes.
   - Run *differentialExpressionAnalysis.R* to identify differentially expressed genes (DEGs).
 
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