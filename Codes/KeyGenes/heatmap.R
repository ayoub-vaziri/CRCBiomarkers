library(data.table)
library(pheatmap)
library(dplyr)

# Set the current working directory to the project path
setwd("PROJECT_PATH")

#### Heat map for validation set ####
#####################################
expr <- as.data.frame(fread("Results/DataProcessing/trainTestSplit/training_set.csv"))
bio <- expr$group

# Gene_sample heatmap
rownames(expr) <- expr$V1
expr <- expr[,-c(1,2)]
expr <- as.data.frame(t(expr))

bio <- gsub(1, "N", bio)
bio <- gsub(2, "C", bio)
bio <- as.factor(bio)
levels(bio) <- c("Tumor", "Normal")

exprData <- expr
bioData <- bio

DEGs <- fread("Results/KeyGenes/differentialExpressionAnalysis/updown.txt", header = FALSE)$V1

data1 <- exprData[DEGs,] %>% na.omit()
data1 <- as.data.frame(t(scale(t(data1))))

samples <- data.frame(Group=as.character(bioData))
rownames(samples) <- colnames(exprData)

n <- rownames(samples)[which(samples$Group == "Normal")]
t <- rownames(samples)[which(samples$Group == "Tumor")]

col_order <- c(n, t)

upgenes <- read.table("Results/KeyGenes/differentialExpressionAnalysis/up.txt")
downgenes <- read.table("Results/KeyGenes/differentialExpressionAnalysis/down.txt")

colnames(upgenes) <- "gene"
colnames(downgenes) <- "gene"

up <- which(rownames(data1) %in% upgenes$gene)
down <- which(rownames(data1) %in% downgenes$gene)

genes <- data.frame(Regulated=rep("X", nrow(data1)))
genes$Regulated[up] <- "up"
genes$Regulated[down] <- "down"
rownames(genes) <- rownames(data1)

p <- pheatmap(data1[, col_order],
              annotation_col = samples,
              annotation_row = genes,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              show_colnames = FALSE,
              show_rownames = FALSE,
              annotation_names_row = FALSE,
              annotation_names_col = FALSE,
              clustering_distance_rows = "correlation",
              clustering_distance_cols = "correlation",
              border_color = NA,
              color = colorRampPalette(c("blue", "white", "red"))(100),
              fontsize = 9,
              fontsize_col = 6
)

png("Results/KeyGenes/heatmap/heatmap.png", height = 2500, width = 3200, res = 300)
p
dev.off()
#####################################