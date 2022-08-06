setwd("C:/Users/ayoub/OneDrive/Desktop/")

interactions <- read.delim("PPI/string_interaction_short_node.csv", sep = ",")

overlap.up <- read.table("DEGs/overlap/overlap.up.txt"); colnames(overlap.up) <- "gene"
overlap.down <- read.table("DEGs/overlap/overlap.down.txt"); colnames(overlap.down) <- "gene"

up <- which(interactions$name %in% overlap.up$gene)
down <- which(interactions$name %in% overlap.down$gene)

interactions$exprType <- "NA"
interactions$exprType[up] <- "up"
interactions$exprType[down] <- "down"

write.table(interactions, file = "PPI/string_interactions_short.node.csv", sep = ",", row.names = F)
