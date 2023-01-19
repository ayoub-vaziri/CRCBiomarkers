library(ggvenn)
library(data.table)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

module <- as.data.frame(fread("Res/CEMiTool/Tables/module.tsv"))

CRGs <- module$genes[which(module$modules %in% c("M1", "M4", "Not.Correlated"))]
DEGs <- fread("Res/limma/updown.txt", header = FALSE)$V1

write.table(CRGs, file = "Res/CEMiTool/CRGs.txt", quote = F, row.names = F, col.names = F)

lstDEGsCRGs <- list(
  DEGs = DEGs,
  CRGs = CRGs
)

png("Res/venn/DEGsCRGs.png", height = 1600, width = 1700, res = 300)
ggvenn(data = lstDEGsCRGs, 
       show_elements = FALSE, 
       show_percentage = FALSE,
       label_sep = "\n", 
       fill_color = c("blue", "red"),
       fill_alpha = 0.3,
       text_size = 3.2,
       set_name_size = 4,
       text_color = "black",
       stroke_alpha = 0.5,
       stroke_size = 0.1,
)
dev.off()

DEGsCRGs <- intersect(DEGs, CRGs)
write.table(DEGsCRGs, file = "Res/PPI/DEGsCRGs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

lassoGenes <- fread("Res/LASSO/lassoGenes.txt", header = FALSE)$V1
rfGenes <- fread("Res/RF/rfGenes.txt", header = FALSE)$V1

lstLassoRf <- list(
  LASSO = lassoGenes,
  `RF-RFE` = rfGenes
)

png("Res/venn/lassoRfGenes.png", height = 1600, width = 1700, res = 300)
ggvenn(data = lstLassoRf, 
       show_elements = TRUE, 
       show_percentage = FALSE,
       label_sep = "\n", 
       fill_color = c("blue", "red"),
       fill_alpha = 0.3,
       text_size = 2.5,
       set_name_size = 3.5,
       text_color = "black",
       stroke_alpha = 0.5,
       stroke_size = 0.1
       )
dev.off()

lassoRFGenes <- intersect(lassoGenes, rfGenes)
write.table(lassoRFGenes, file = "Res/venn/lassoRfGenes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
