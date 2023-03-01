library(pROC)
library(data.table)

set.seed(123)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### ROC curve function ####
############################
rocCurve <- function(data, genes, path, mfrow) {
  png(filename = path, width = 3600, height = 1800, res = 300)
  par(mfrow=mfrow)
  for(g in genes) {
    ROC <- roc(as.formula(paste("status ~", g)), data = data, auc = TRUE, ci = TRUE)
    plot.roc(ROC, 
             print.auc = T, 
             col = "red",
             ci = T,
             print.auc.pattern = "       AUC: %.3f \n95%% CI: %.3f-%.3f",
             print.auc.x = 0.64,
             print.auc.y = 0.49,
             legacy.axes = TRUE,
             main = g,
             cex.main = 1.6,
             lwd = 1.1,
             print.auc.cex = 1.3,
             cex.lab = 1.3
    )
  }
  dev.off()
}
############################

#### ROC curve for validation data ####
#######################################
expr_data <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
bio_data <- fread("Res/exprBatchBioData/bioColorectal.txt", header = FALSE)$V1

lassoRfGenes <- fread("Res/venn/lassoRfGenes.txt", header = FALSE)$V1
lassoRfGenes <- sort(lassoRfGenes)

valid_smps <- 1:102

valid_expr <- expr_data[valid_smps,]
rownames(valid_expr) <- valid_expr$V1
valid_expr <- valid_expr[,-1]
valid_expr <- valid_expr[,lassoRfGenes]

valid_bio <- bio_data[valid_smps]
valid_bio <- as.factor(valid_bio)
levels(valid_bio) <- c("Normal", "Tumor")

valid_data <- cbind(status=valid_bio, valid_expr)

rocCurve(valid_data, lassoRfGenes, "Res/supplementary/ROC_valid.png", c(2, 5))
rocCurve(valid_data, lassoRfGenes[-7], "Res/ROC/ROC_valid.png", c(2, 4))
#######################################

#### ROC curve for train data ####
##################################
train_smps <- 103:602

train_expr <- expr_data[train_smps,]
rownames(train_expr) <- train_expr$V1
train_expr <- train_expr[,-1]
train_expr <- train_expr[,lassoRfGenes]

train_bio <- bio_data[train_smps]
train_bio <- as.factor(train_bio)
levels(train_bio) <- c("Normal", "Tumor")

train_data <- cbind(status=train_bio, train_expr)

rocCurve(train_data, lassoRfGenes, "Res/supplementary/ROC_train.png", c(2, 5))
rocCurve(train_data, lassoRfGenes[-7], "Res/ROC/ROC_train.png", c(2, 4))
##################################