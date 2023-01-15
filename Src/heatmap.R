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

smps <- 1:52

exprValid <- expr[,smps]
bioValid <- bio[smps]

DEGsRGs <- fread("Res/PPI/DEGsRGs.txt", header = FALSE)$V1

data1 <- exprValid[DEGsRGs,] %>% na.omit()
data1 <- as.data.frame(t(scale(t(data1))))

samples <- data.frame(Group=as.character(bioValid))
rownames(samples) <- colnames(exprValid)

n <- rownames(samples)[which(samples$Group == "Normal")]
t <- rownames(samples)[which(samples$Group == "Tumor")]

col_order <- c(n, t)

upgenes <- read.table("Res/limma/up.txt")
downgenes <- read.table("Res/limma/down.txt")

colnames(upgenes) <- "gene"
colnames(downgenes) <- "gene"

up <- which(rownames(data1) %in% upgenes$gene)
down <- which(rownames(data1) %in% downgenes$gene)

genes <- data.frame(`Regulated_genes`=rep("X", nrow(data1)))
genes$Regulated_genes[up] <- "up"
genes$Regulated_genes[down] <- "down"
rownames(genes) <- rownames(data1)

p <- pheatmap(data1[, col_order],
              annotation_col = samples,
              annotation_row = genes,
              cluster_cols = TRUE,
              cluster_rows = TRUE,
              show_colnames = TRUE,
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

png("Res/heatmap/heatmap_validation_set.png", height = 2000, width = 2400, res = 300)
p
dev.off()


# Samples heatmap
data2 <- exprValid[DEGsRGs,] %>% na.omit()
data2 <- as.data.frame(t(exprValid))
data2 <- as.matrix(scale(data2, center = FALSE))
data2 <- cbind(STATUS=as.character(bioValid), data2)
  
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
         border_color = NA,
         fontsize = 9,
         fontsize_row = 6,
         fontsize_col = 6
         
)
dev.off()
#####################################