library(MLmetrics)
library(pROC)
library(randomForest)
library(caret)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### ML models evaluation ####
##############################
# Load trained RF and SVM models
RF <- readRDS("Res/trainEvaluateML/RF.rds")
SVM <- readRDS("Res/trainEvaluateML/SVM.rds")

# log2 transform
log2trans <- function(expr) {
  quan <- as.numeric(quantile(expr, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (quan[5] > 100) || (quan[6]-quan[1] > 50 && qaun[2] > 0)
  if(LogC) { 
    # expr[which(expr <= 0)] <- NaN
    expr <- log2(expr+1) 
  }
  return(expr)
}

pred <- function(exprData, model) {
  x <- log2trans(exprData)
  y.pred <- as.character(predict(model, x))
  y.prob <- as.data.frame(predict(model, x, type = "prob"))
  return(list(y.pred, y.prob))
}

rf.metrics <- function(x, y, accession, model) {
  y.true <- as.factor(y)
  y.pred <- as.factor(pred(x, model)[[1]])
  y.prob <- as.data.frame(pred(x, model)[[2]])
  roc <- roc(y.true ~ as.numeric(y.prob$Tumor))
  cm <- confusionMatrix(y.true, y.pred, positive = "Tumor")$table
  write.csv(cm, paste0("Res/externalTesting/", accession, ".ConfusionMatrix.csv"))
  mats <- data.frame("AUC"=roc$auc[1],
                     "accuracy"=Accuracy(y.pred, y.true),
                     "sensitivity"=Sensitivity(y.true, y.pred, positive = "Tumor"),
                     "specificity"=Specificity(y.true, y.pred, positive = "Tumor"),
                     "precision"=Precision(y.true, y.pred, positive = "Tumor"),
                     "recall"=Recall(y.true, y.pred, positive = "Tumor"),
                     "f1score"=F1_Score(y.true, y.pred, positive = "Tumor")
                     )

  write.csv(mats, paste0("Res/externalTesting/", accession, ".Metrics.csv"))
  return(mats)
}


svm.metrics <- function(x, y, accession, model) {
  y.true <- as.factor(y)
  y.pred <- unname(predict(model, newdata = x))
  y.prob <- as.data.frame(attr(predict(model, newdata = x, probability = TRUE), "probabilities"))
  
  roc <- roc(y.true ~ as.numeric(y.prob$Tumor))
  
  cm <- table(y.true, y.pred)
  write.csv(cm, paste0("Res/externalTesting/", accession, ".ConfusionMatrix.csv"))
  
  mats <- data.frame(
                    "AUC"=roc$auc[1],
                    "accuracy"=Accuracy(y.pred, y.true),
                    "sensitivity"=Sensitivity(y.true, y.pred, positive = "Tumor"),
                    "specificity"=Specificity(y.true, y.pred, positive = "Tumor"),
                    "precision"=Precision(y.true, y.pred, positive = "Tumor"),
                    "recall"=Recall(y.true, y.pred, positive = "Tumor"),
                    "f1score"=F1_Score(y.true, y.pred, positive = "Tumor")
                    )
  write.csv(mats, paste0("Res/externalTesting/", accession, ".Metrics.csv"))
  return(mats)
}
##############################


#### Evaluation of ML models on external testing set ####
#########################################################
gse21815 <- as.data.frame(fread("Res/externalTesting/GSE21815.csv"))
rownames(gse21815) <- gse21815$V1
gse21815 <- gse21815[,-1]
gse21815[,-1] <- log2(gse21815[,-1]+1)

gse106582 <- as.data.frame(fread("Res/externalTesting/GSE106582.csv"))
rownames(gse106582) <- gse106582$V1
gse106582 <- gse106582[,-1]

mergedData <- rbind(gse21815, gse106582)

x <- mergedData[,-1]
y <- mergedData[,1]

rf.metrics(x, y, "External testing set for RF", RF)
svm.metrics(x, y, "External testing set for SVM", SVM)
#########################################################