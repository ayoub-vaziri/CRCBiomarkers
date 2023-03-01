library(data.table)
library(factoextra)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### PCA plot for validation set ####
#####################################
train_data <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
bio_data <- fread("Res/exprBatchBioData/bioColorectal.txt", header = FALSE)$V1

DEGsCRG <- fread("Res/PPI/DEGsCRGs.txt", header = FALSE)$V1

expr_data <- train_data[,DEGsCRG]
rownames(expr_data) <- expr_data$V1
expr_data <- expr_data[,-1]

bio_data <- as.factor(bio_data)
levels(bio_data) <- c("Normal", "Tumor")
bio_data <- as.character(bio_data)

expr_data <- cbind(STATUS=bio_data, expr_data)

pca <- prcomp(expr_data[,-1], scale. = FALSE, center = TRUE)

ppca <- fviz_pca_ind(pca,
                  label = "none",
                  palette = c("cyan3", "coral2"),
                  habillage = expr_data$STATUS,
                  addEllipses = TRUE,
                  ellipse.level=0.75)
  
png("Res/PCA/pcaplot.png", height = 2000, width = 2400, res = 300)
ppca +
  theme_bw() +
  labs(title = "", x = "PC1", y = "PC2")
dev.off()
#####################################