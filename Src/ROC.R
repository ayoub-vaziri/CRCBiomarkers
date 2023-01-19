library(dplyr)
library(ggpubr)
library(data.table)
library(pROC)

set.seed(123)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### ROC curve function ####
############################
rocCurve <- function(data, genes, path) {
  png(filename = path, width = 3600, height = 1800, res = 300)
  par(mfrow=c(2, 4))
  for(g in genes) {
    ROC <- roc(as.formula(paste("status ~", g)), data = data, auc = TRUE, ci = TRUE)
    plot.roc(ROC, print.auc = T, 
             col = "red",
             ci = T,
             print.auc.pattern = "       AUC: %.3f \n95%% CI: %.3f-%.3f",
             print.auc.x = 0.64,
             print.auc.y = 0.49,
             legacy.axes = TRUE,
             main = g,
             lwd = 1.1,
             print.auc.cex = 1.2,
             cex.lab = 1.3
    )
  }
  dev.off()
}
############################

#### ROC curve for validation data ####
#######################################
train_data <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
bio_data <- fread("Res/exprBatchBioData/bioColorectal.txt", header = FALSE)$V1
lassoRfGenes <- fread("Res/venn/lassoRfGenes.txt", header = FALSE)$V1

smps <- 1:52

valid_data <- train_data[smps,]
rownames(valid_data) <- valid_data$V1
valid_data <- valid_data[,-1]
valid_data <- valid_data[,lassoRfGenes]

bio_data <- bio_data[smps]
bio_data <- as.factor(bio_data)
levels(bio_data) <- c("Normal", "Tumor")

valid_data <- cbind(status=bio_data, valid_data)

rocCurve(valid_data, lassoRfGenes, "Res/supplementary/ROC_valid.png")
rocCurve(valid_data, lassoRfGenes[-1], "Res/ROC/ROC_valid.png")
#######################################

#### ROC curve for train data ####
##################################
rownames(train_data) <- train_data$V1
train_data <- train_data[,-1]

status <- read.table("Res/exprBatchBioData/bioColorectal.txt", header = FALSE)$V1
status <- as.factor(status)
levels(status) <- c(0, 1)
status <- as.numeric(as.character(status))

train_data <- cbind(STATUS=status, train_data)

rocCurve(train_data, lassoRfGenes, "Res/supplementary/ROC_train.png")
rocCurve(train_data, lassoRfGenes[-1], "Res/ROC/ROC_train.png")
##################################