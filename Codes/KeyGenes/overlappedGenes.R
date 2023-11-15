library(ggvenn)
library(data.table)

# Set the current working directory to the project path
setwd("PROJECT_PATH")

degs <- fread("Results/KeyGenes/differentialExpressionAnalysis/updown.txt", header = FALSE)$V1
modules <- fread("Results/KeyGenes/geneCoexpressionAnalysis/modulesGenes.txt", header = FALSE)$V1

lstDegsModules <- list(
  DEGs = degs,
  M1 = modules
)

png("Results/KeyGenes/overlappedGenes/overlappedGenes.png", height = 1600, width = 1600, res = 300)
ggvenn(data = lstDegsModules, 
       show_elements = FALSE, 
       show_percentage = FALSE,
       label_sep = "\n", 
       fill_color = c("blue", "red"),
       fill_alpha = 0.3,
       text_size = 6,
       set_name_size = 5.5,
       text_color = "black",
       stroke_alpha = 0.5,
       stroke_size = 0.2,
)
dev.off()

overlappedGenes <- intersect(degs, modules)
write.table(overlappedGenes, file = "Results/KeyGenes/overlappedGenes/overlappedGenes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)