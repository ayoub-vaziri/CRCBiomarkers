library(data.table)
library(dplyr)
library(factoextra)
library(pheatmap)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### PCA plot for validation set ####
#####################################
train_data <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
bio_data <- fread("Res/exprBatchBioData/bioColorectal.txt", header = FALSE)$V1

smps <- 1:52

DEGsCRG <- fread("Res/PPI/DEGsCRGs.txt", header = FALSE)$V1

valid_data <- train_data[smps,DEGsCRG]
rownames(valid_data) <- valid_data$V1
valid_data <- valid_data[,-1]

bio_data <- bio_data[smps]
bio_data <- as.factor(bio_data)
levels(bio_data) <- c("Normal", "Tumor")
bio_data <- as.character(bio_data)

valid_data <- cbind(STATUS=bio_data, valid_data)

pca <- prcomp(valid_data[,-1], scale. = T)

ppca <- fviz_pca_ind(pca,
                  label = "none",
                  palette = c("cyan3", "coral2"),
                  habillage = valid_data$STATUS,
                  addEllipses = TRUE,
                  ellipse.level=0.8)
  
png("Res/PCA/pcaplot.png", height = 1600, width = 2000, res = 300)
ppca +
  theme_bw() +
  labs(title = "", x = "PC1", y = "PC2")
dev.off()
#####################################