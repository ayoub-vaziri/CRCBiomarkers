library(data.table)
library(pheatmap)
library(dplyr)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Heat map for validation set ####
#####################################
expr <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
bio <- fread("Res/exprBatchBioData/bioColorectal.txt")$V1

# Gene_sample heatmap
rownames(expr) <- expr$V1
expr <- expr[,-1]
expr <- as.data.frame(t(expr))

bio <- gsub(1, "N", bio)
bio <- gsub(2, "C", bio)
bio <- as.factor(bio)
levels(bio) <- c("Tumor", "Normal")

exprData <- expr
bioData <- bio

DEGsRGs <- fread("Res/PPI/DEGsCRGs.txt", header = FALSE)$V1

data1 <- exprData[DEGsRGs,] %>% na.omit()
data1 <- as.data.frame(t(scale(t(data1))))

samples <- data.frame(Group=as.character(bioData))
rownames(samples) <- colnames(exprData)

n <- rownames(samples)[which(samples$Group == "Normal")]
t <- rownames(samples)[which(samples$Group == "Tumor")]

col_order <- c(n, t)

upgenes <- read.table("Res/limma/up.txt")
downgenes <- read.table("Res/limma/down.txt")

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

png("Res/heatmap/heatmap.png", height = 2500, width = 3200, res = 300)
p
dev.off()

# Top20 heatmap
top20 <- fread("Res/PPI/top20.txt", header = FALSE)$V1

dataTop20 <- exprData[top20,] %>% na.omit()
dataTop20 <- as.data.frame(t(scale(t(dataTop20))))

up <- which(rownames(dataTop20) %in% upgenes$gene)
down <- which(rownames(dataTop20) %in% downgenes$gene)

genes <- data.frame(Regulated=rep("X", nrow(dataTop20)))
genes$Regulated[up] <- "up"
genes$Regulated[down] <- "down"
rownames(genes) <- rownames(dataTop20)

p <- pheatmap(dataTop20[, col_order],
              annotation_col = samples,
              annotation_row = genes,
              cluster_cols = FALSE,
              cluster_rows = FALSE,
              show_colnames = FALSE,
              show_rownames = TRUE,
              annotation_names_row = FALSE,
              annotation_names_col = FALSE,
              clustering_distance_rows = "correlation",
              clustering_distance_cols = "correlation",
              border_color = NA,
              color = colorRampPalette(c("blue", "white", "red"))(100),
              fontsize = 9,
              fontsize_row = 13
)

png("Res/heatmap/heatmap_top20.png", height = 1400, width = 3400, res = 300)
p
dev.off()

# Samples heatmap
data2 <- exprData[DEGsRGs,] %>% na.omit()
data2 <- as.data.frame(t(exprData))
data2 <- as.matrix(scale(data2, center = FALSE))
data2 <- cbind(STATUS=as.character(bioData), data2)
  
samples <- data.frame(Group=data2[,1])
rownames(samples) <- rownames(data2)

data2 <- as.data.frame(t(data2[,-1]))

for(x in 1:ncol(data2)) data2[,x] <- as.numeric(data2[,x])

png("Res/heatmap/heatmapSamples.png", height = 2000, width = 2400, res = 300)
pheatmap(cor(data2),
         annotation_col = samples,
         annotation_row = samples,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         cutree_rows = TRUE,
         cutree_cols = TRUE,
         show_colnames = FALSE,
         show_rownames = FALSE,
         #border_color = NA,
         fontsize = 9,
         fontsize_row = 6,
         fontsize_col = 6
)
dev.off()
#####################################