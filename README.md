# CRCBiomarkers

## Overview
This github repository contains the data files and analysis code used for the scientific paper titled **"Identifying novel gene biomarkers for colorectal cancer diagnosis through bioinformatics analysis and machine learning"**.
The files are organised into four folders:
* _Data_: which contains all the transcriptomic data required to perform the analyses described in the paper
* _Codes_: contains the R code to reproduce all  analyses.
* _Results_: contains all the results produced by the R scripts.
* _CRCTest_: contains a trained random forest model and code for diagnosing tumor samples from normal samples.

## Reproducing the results
This repository contains all the code necessary to reproduce the results in the paper:
 - **Codes/DataProcessing** contains the code necessary to pre-process dataset.
   - Run *mergeDatasets.R* to merge the five microarray datasets (GSE10950, GSE25070, GSE41328, GSE74602, and GSE142279) based on their common genes.

## Required software
The scripts use core R functionality and several publicly available R packages listed below. Version numbers in brackets correspond to the versions of the packages that were used to develop and debug these scripts.

 - R (4.3.1)
 - RStudio (2022.07.2+576): Optional, testing functions and running the code step-by-step.
 - 