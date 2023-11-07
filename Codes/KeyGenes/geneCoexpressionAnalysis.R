library(data.table)
library(CEMiTool)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load and process data ####
###############################
trainingSet <- as.data.frame(fread("Res/trainTestSplit/training_set.csv"))

exprData <- trainingSet
phenData <- trainingSet$group

gmtData <- read_gmt(system.file("extdata", "pathways.gmt", package = "CEMiTool"))
intData <- read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))

rownames(exprData) <- exprData$V1
exprData <- exprData[,-c(1,2)]
exprData <- as.data.frame(t(exprData))

write.csv(exprData, file = paste0("Res/geneCoexpressionAnalysis/GeneExpressionProfileOf",ncol(exprData),"Samples.csv"), quote = F, row.names = T)

phenData <- data.frame(ID=colnames(exprData), Group=phenData)
phenData$Group <- gsub(1, "Normal", phenData$Group)
phenData$Group <- gsub(2, "Tumor", phenData$Group)

write.csv(phenData, file = paste0("Res/geneCoexpressionAnalysis/DiseasePhenotypesOf",ncol(exprData),"Samples.csv"), quote = F, row.names = T)
###############################

#### WGCNA using CEMiTool package ####
######################################
cem <- cemitool(expr = exprData, 
                annot = phenData,
                gmt = gmtData,
                interactions = intData,
                sample_name_column = "ID", 
                class_column = "Group",
                plot = TRUE,
                verbose = TRUE
                )

generate_report(cem, directory = "Res/geneCoexpressionAnalysis/Report", force = T, output_format = "html_document")
write_files(cem, directory = "Res/geneCoexpressionAnalysis/Tables", force = T)
save_plots(cem, value = c("all"), force = T, directory = "Res/geneCoexpressionAnalysis/Plots")



png("Res/geneCoexpressionAnalysis/beta_r2.png", height = 1600, width = 1600, res = 300)
show_plot(cem, "beta_r2")
dev.off()

png("Res/geneCoexpressionAnalysis/mean_k.png", height = 1600, width = 1600, res = 300)
show_plot(cem, "mean_k")
dev.off()

png("Res/geneCoexpressionAnalysis/gsea.png", height = 1600, width = 1600, res = 300)
show_plot(cem, "gsea")
dev.off()

png("Res/geneCoexpressionAnalysis/ora_M1.png", height = 1600, width = 1600, res = 300)
show_plot(cem, "ora")$M1
dev.off()

mods <- fread("Res/geneCoexpressionAnalysis/Tables/module.tsv")
sigMods <- c("M1")
sigGenes <- mods$genes[which(mods$modules %in% sigMods)]
write.table(sigGenes, "Res/geneCoexpressionAnalysis/modulesGenes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
