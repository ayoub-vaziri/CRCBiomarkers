# load packages
library(CEMiTool)
library(magrittr)
library(doParallel)

# set working directory
setwd("C:/Users/ayoub/OneDrive/Desktop/CEMiTool/TCGA_COADREAD/")

#### load data ####
clinData <- read.table("Data/COADREAD.clin.merged.picked.txt", sep = "\t", header = T)
exprData <- read.table("Data/COADREAD.uncv2.mRNAseq_RSEM_normalized_log2.txt", sep = "\t", header = T)

#### data preprocessing ####
exprData$gene <- exprData$gene %>% gsub("\\|.*", "", .)
exprData <- exprData[-which(exprData$gene == "?"),] # removing no symbols

n_occur <- data.frame(table(exprData$gene))
dup_id <- which(exprData$gene %in% n_occur$Var1[n_occur$Freq > 1])
exprData <- exprData[-dup_id[2],] # removing duplicate symbol

exprData <- exprData[,c(1,which(substr(colnames(exprData), 14, 15) == "01"))]
gene_symbol <- exprData$gene

colnames(clinData)[2:length(colnames(clinData))] <- paste0(toupper(colnames(clinData)), ".01")[2:length(colnames(clinData))]

int.samples <- Reduce(intersect, list(colnames(clinData), colnames(exprData)))

exprData <- exprData[,which(colnames(exprData) %in% int.samples)]
rownames(exprData) <- gene_symbol
exprData <- na.omit(exprData)

clin.ids <- unlist(lapply(colnames(exprData), function(x) which(colnames(clinData) == x)))
clinData <- clinData[3,clin.ids]
clinData <- as.data.frame(t(clinData))
clinData <- data.frame(sample_id=rownames(clinData), vital_status=clinData$`3`) %>% na.omit()
rownames(clinData) <- NULL
clinData$vital_status <- gsub("1", "Dead", clinData$vital_status)
clinData$vital_status <- gsub("0", "Live", clinData$vital_status)

exprData <- exprData[, which(colnames(exprData) %in% clinData$sample_id)]

#### WGCNA using CEMiTool package ####
expr0 <- exprData
sample_annot <- clinData

registerDoParallel()

cem <- cemitool(expr = expr0, 
                annot = sample_annot, 
                sample_name_column = "sample_id", 
                class_column = "vital_status", 
                verbose = T, 
                plot = T)

generate_report(cem, directory = "Res/Report", force = T, output_format = "html_document")
write_files(cem, directory = "Res/Tables", force = T)
save_plots(cem, value = c("all"), force = T, directory = "Res/Plots")
