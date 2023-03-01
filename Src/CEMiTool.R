library(data.table)
library(CEMiTool)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load and process data ####
###############################
exprData <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
phenData <- fread("Res/exprBatchBioData/bioColorectal.txt")$V1

gmtData <- read_gmt(system.file("extdata", "pathways.gmt", package = "CEMiTool"))
intData <- read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))

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
                gmt = gmtData,
                interactions = intData,
                sample_name_column = "ID", 
                class_column = "STATUS",
                merge_similar = TRUE,
                plot = TRUE,
                verbose = TRUE)

generate_report(cem, directory = "Res/CEMiTool/Report", force = T, output_format = "html_document")
write_files(cem, directory = "Res/CEMiTool/Tables", force = T)
save_plots(cem, value = c("all"), force = T, directory = "Res/CEMiTool/Plots")

png("Res/CEMiTool/beta_r2.png", height = 1600, width = 1600, res = 300)
show_plot(cem, "beta_r2")
dev.off()

png("Res/CEMiTool/mean_k.png", height = 1600, width = 1600, res = 300)
show_plot(cem, "mean_k")
dev.off()

png("Res/CEMiTool/gsea.png", height = 1600, width = 1600, res = 300)
show_plot(cem, "gsea")
dev.off()

options(ggrepel.max.overlaps = 18)

png("Res/CEMiTool/M1_interaction.png", height = 1600, width = 1600, res = 300)
plts <- plot_interactions(cem)
plots <- show_plot(plts, "interaction")
plots[1]
dev.off()

png("Res/CEMiTool/M44_interaction.png", height = 1600, width = 1600, res = 300)
plts <- plot_interactions(cem, n=15)
plots <- show_plot(plts, "interaction")
plots[4]
dev.off()
######################################