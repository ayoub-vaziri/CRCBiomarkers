library(dplyr)
library(ggpubr)
library(data.table)
library(pROC)

set.seed(123)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### ROC curve for validation data ####
#######################################
train_data <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
bio_data <- fread("Res/exprBatchBioData/bioColorectal.txt", header = FALSE)$V1
lassoGenes <- fread("Res/LASSO/lassoGenes.txt", header = FALSE)$V1

smps <- 1:52

valid_data <- train_data[smps,]
rownames(valid_data) <- valid_data$V1
valid_data <- valid_data[,-1]
valid_data <- valid_data[,lassoGenes]

bio_data <- bio_data[smps]
bio_data <- as.factor(bio_data)
levels(bio_data) <- c("Normal", "Tumor")

valid_data <- cbind(status=bio_data, valid_data)

# ROC curve
png(filename = "Res/ROC/ROC_valid.png", width = 2800, height = 2800, res = 300)
par(mfrow=c(3, 3))
for(g in lassoGenes) {
  ROC <- roc(as.formula(paste("status ~", g)), data = valid_data, auc = TRUE, ci = TRUE)
  plot.roc(ROC, print.auc = T, 
           col = "red",
           ci = T,
           print.auc.pattern = "        AUC: %.3f \n95%% CI: %.3f-%.3f",
           print.auc.x = 0.57,
           print.auc.y = 0.49,
           legacy.axes = TRUE,
           main = g,
           lwd = 1.1,
           print.auc.cex = 1.2,
           cex.lab = 1.3
           
  )
}
dev.off()
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

# ROC curve
png(filename = "Res/ROC/ROC_train.png", width = 2800, height = 2800, res = 300)
par(mfrow=c(3, 3))
for(g in lassoGenes) {
  ROC <- roc(as.formula(paste("STATUS ~", g)), data = train_data, auc = TRUE, ci = TRUE)
  plot.roc(ROC, print.auc = T, 
           col = "red",
           ci = T,
           print.auc.pattern = "        AUC: %.3f \n95%% CI: %.3f-%.3f",
           print.auc.x = 0.57,
           print.auc.y = 0.49,
           legacy.axes = TRUE,
           main = g,
           lwd = 1.1,
           print.auc.cex = 1.2,
           cex.lab = 1.3
  )
}
dev.off()
##################################