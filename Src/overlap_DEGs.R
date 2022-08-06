library(data.table)

setwd("C:/Users/ayoub/OneDrive/Desktop/")

data <- c("GSE39582", "GSE110223", "GSE110224", "TCGA_COADREAD")

down.list <- list()
up.list <- list()

for(d in data) {
  down.list[[d]] <- fread(paste0("DEGs/up.down/", d, ".downregulated.txt"), header = F)$V1
  up.list[[d]] <- fread(paste0("DEGs/up.down/", d, ".upregulated.txt"), header = F)$V1
}

overlap.down <- Reduce(intersect, down.list)
overlap.up <- Reduce(intersect, up.list)
overlap.updown <- union(overlap.down, overlap.up)
  
write.table(overlap.down, file = "DEGs/overlap/overlap.down.txt", quote = F, row.names = F, col.names = F)
write.table(overlap.up, file = "DEGs/overlap/overlap.up.txt", quote = F, row.names = F, col.names = F) 
write.table(overlap.updown, file = "DEGs/overlap/overlap.updown.txt", quote = F, row.names = F, col.names = F) 
