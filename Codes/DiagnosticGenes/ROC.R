library(pROC)
library(data.table)

set.seed(123)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### ROC curve function ####
############################
rocCurve <- function(data, genes, path, mfrow, auc.x, auc.y, w, h) {
  png(filename = path, width = w, height = h, res = 300)
  par(mfrow=mfrow)
  for(g in genes) {
    roc.train <- roc(as.formula(paste("group ~", g)), data = data[[1]], auc = TRUE, ci = TRUE)
    roc.test <- roc(as.formula(paste("group ~", g)), data = data[[2]], auc = TRUE, ci = TRUE)
    plot.roc(roc.train,
             col = "#1F77B4",
             legacy.axes = TRUE,
             main = g,
             cex.main = 1.2,
             lwd = 2,
             cex.lab = 1.1
    )
    lines(roc.test, col="salmon")
    legend("bottomright",
           legend=c(paste("train AUC:", signif(roc.train$auc*100, digits = 4), "%"),
                    paste("test AUC: ", signif(roc.test$auc*100, digits = 4), "%")),
           col=c("#1F77B4", "salmon"),
           lty = 1,
           lwd = 2
           )
  }
  dev.off()
}
############################


lassoGenes <- fread("Res/LASSO/lassoGenes.txt", header = FALSE)$V1
lassoGenes <- sort(lassoGenes)


#### ROC curve for train data ####
##################################
trainingSet <- as.data.frame(fread("Res/trainTestSplit/training_set.csv"))

rownames(trainingSet) <- trainingSet$V1
group <- trainingSet$group
group <- as.factor(group)
levels(group) <- c("Normal", "Tumor")
trainingSet <- trainingSet[,-c(1,2)]

trainingSet <- trainingSet[, lassoGenes]
train <- cbind(group=group, trainingSet)
##################################


#### ROC curve for testing data ####
####################################
testSet <- as.data.frame(fread("Res/trainTestSplit/testing_set.csv"))

rownames(testSet) <- testSet$V1

group <- testSet$group
group <- as.factor(group)
levels(group) <- c("Normal", "Tumor")
testSet <- testSet[,-c(1,2)]

testSet <- testSet[, lassoGenes]
test <- cbind(group=group, testSet)
####################################


#### ROC curve visualization ####
#################################
trainTestData <- list(train, test)
rocCurve(trainTestData, lassoGenes, "Res/supplementary/ROC.png", c(2, 5), 0.78, 0.50, 3600, 1450)
rocCurve(trainTestData, lassoGenes[-c(3,10)], "Res/ROC/ROC.png", c(2, 4), 0.78, 0.50, 3000, 1500)
#################################