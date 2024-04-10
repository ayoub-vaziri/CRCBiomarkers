library(ggvenn)
library(data.table)

# Set the current working directory to the project path
setwd(project_path)

degs <- fread("Results/differentialExpressionAnalysis/updown.txt", header = FALSE)$V1
modules <- fread("Results/geneCoexpressionAnalysis/modulesGenes.txt", header = FALSE)$V1

lstDegsModules <- list(
  DEGs = degs,
  M1 = modules
)

png("Results/overlappedGenes/DEGsModules.png", height = 1600, width = 1600, res = 300)
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

DEGsModules <- intersect(degs, modules)
write.table(DEGsModules, file = "Results/overlappedGenes/DEGsModules.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
