library(data.table)

top.tables <- c("GSE39582", "GSE110223", "GSE110224", "TCGA_COADREAD")

setwd("C:/Users/ayoub/OneDrive/Desktop/DEGs/")

for(tt in top.tables) {
  tT <- as.data.frame(fread(paste0("top.table/", tt, ".top.table.tsv")))
  
  up.regulated <- subset(tT, tT$adj.P.Val < 0.05 & tT$logFC > 1)
  up.regulated <- unique(up.regulated$Gene.symbol)
  up.regulated <- strsplit(up.regulated, split = '///')
  up.regulated <- unique(unlist(up.regulated))
  write.table(up.regulated, file = paste0("up.down/", tt, ".upregulated.txt"), row.names = F, col.names = F, quote = F)
  
  down.regulated <- subset(tT, tT$adj.P.Val < 0.05 & tT$logFC < -1)
  down.regulated <- unique(down.regulated$Gene.symbol)
  down.regulated <- strsplit(down.regulated, split = '///')
  down.regulated <- unique(unlist(down.regulated))
  write.table(down.regulated, file = paste0("up.down/", tt, ".downregulated.txt"), row.names = F, col.names = F, quote = F)
}
