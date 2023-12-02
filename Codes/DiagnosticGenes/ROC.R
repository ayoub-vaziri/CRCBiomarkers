library(pROC)
library(data.table)

set.seed(123)

# Set the current working directory to the project path
setwd("PROJECT_PATH")

#### ROC curve function ####
############################
rocCurve <- function(data, genes, path, mfrow, auc.x, auc.y, w, h) {
  png(filename = path, width = w, height = h, res = 300)
  par(mfrow=mfrow)
  for(g in genes) {
    roc.train <- roc(as.formula(paste("group ~", g)), data = data[[1]], auc = TRUE, ci = TRUE)
    roc.valid <- roc(as.formula(paste("group ~", g)), data = data[[2]], auc = TRUE, ci = TRUE)
    plot.roc(roc.train,
             col = "#1F77B4",
             legacy.axes = TRUE,
             main = g,
             cex.main = 1.2,
             lwd = 2,
             cex.lab = 1.1
    )
    lines(roc.valid, col="salmon")
    legend("bottomright",
           legend=c(paste("train AUC:", signif(roc.train$auc*100, digits = 4), "%"),
                    paste("valid AUC: ", signif(roc.valid$auc*100, digits = 4), "%")),
           col=c("#1F77B4", "salmon"),
           lty = 1,
           lwd = 2
           )
  }
  dev.off()
}
############################


lassoGenes <- fread("Results/DiagnosticGenes/LASSO/lassoGenes.txt", header = FALSE)$V1
lassoGenes <- sort(lassoGenes)


#### ROC curve for train data ####
##################################
trainingSet <- as.data.frame(fread("Results/DataProcessing/trainvalidSplit/training_set.csv"))

rownames(trainingSet) <- trainingSet$V1
group <- trainingSet$group
group <- as.factor(group)
levels(group) <- c("Normal", "Tumor")
trainingSet <- trainingSet[,-c(1,2)]

trainingSet <- trainingSet[, lassoGenes]
train <- cbind(group=group, trainingSet)
##################################


#### ROC curve for validing data ####
####################################
validSet <- as.data.frame(fread("Results/DataProcessing/trainvalidSplit/validation_set.csv"))

rownames(validSet) <- validSet$V1

group <- validSet$group
group <- as.factor(group)
levels(group) <- c("Normal", "Tumor")
validSet <- validSet[,-c(1,2)]

validSet <- validSet[, lassoGenes]
valid <- cbind(group=group, validSet)
####################################


#### ROC curve visualization ####
#################################
trainValidData <- list(train, valid)
#rocCurve(trainValidData, lassoGenes, "Res/supplementary/ROC.png", c(2, 5), 0.78, 0.50, 3600, 1450)
rocCurve(trainValidData, lassoGenes[-c(3,10)], "Results/DiagnosticGenes/ROC/ROC.png", c(2, 4), 0.78, 0.50, 3000, 1500)
#################################