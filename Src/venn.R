library(VennDiagram)
library(data.table)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

module <- as.data.frame(fread("Res/CEMiTool/Tables/module.tsv"))
RGs <- module$genes[which(module$modules %in% c("M1", "M4", "Not.Correlated"))]

write.table(RGs, file = "Res/CEMiTool/RGs.txt", quote = F, row.names = F, col.names = F)

display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

DEGs <- fread("Res/limma/updown.txt", header = FALSE)$V1

lstDEGsRGs <- list(
  DEGs = DEGs,
  RGs = RGs
)

png("Res/venn/DEGsRGs.png", height = 1600, width = 1700, res = 300)
display_venn(
  lstDEGsRGs,
  category.names = c("DEGs", "RGs"),
  fill = c(DEGs="blue", RGs="red"),
  cex = 1.5,
  cat.cex = 1.2
)
dev.off()

DEGsRGs <- intersect(DEGs, RGs)

write.table(DEGsRGs, file = "Res/PPI/DEGsRGs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
