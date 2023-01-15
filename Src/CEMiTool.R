library(data.table)
library(dplyr)
library(CEMiTool)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load and process data ####
###############################
exprData <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
phenData <- fread("Res/exprBatchBioData/bioColorectal.txt")$V1

rownames(exprData) <- exprData$V1
exprData <- exprData[,-1]
exprData <- as.data.frame(t(exprData))

write.csv(exprData, file = paste0("Res/CEMiTool/GeneExpressionProfileOf",ncol(exprData),"Samples.csv"), quote = F, row.names = T)

phenData <- data.frame(ID=colnames(exprData), STATUS=phenData)
phenData$STATUS <- gsub(1, "Normal", phenData$STATUS)
phenData$STATUS <- gsub(2, "Tumor", phenData$STATUS)

write.csv(phenData, file = paste0("Res/CEMiTool/DiseasePhenotypesOf",ncol(exprData),"Samples.csv"), quote = F, row.names = T)
###############################

#### WGCNA using CEMiTool package ####
######################################
cem <- cemitool(expr = exprData, 
                annot = phenData, 
                sample_name_column = "ID", 
                class_column = "STATUS", 
                merge_similar = T,
                verbose = T, 
                plot = T)

generate_report(cem, directory = "Res/CEMiTool/Report", force = T, output_format = "html_document")
write_files(cem, directory = "Res/CEMiTool/Tables", force = T)
save_plots(cem, value = c("all"), force = T, directory = "Res/CEMiTool/Plots")
######################################