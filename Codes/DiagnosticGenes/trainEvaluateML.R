library(data.table)
library(ggplot2)
library(caret)
library(pROC)
library(MLmetrics)
library(randomForest)
library(e1071)
library(dplyr)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

set.seed(123)

#### Load and process train data ####
#####################################
trainSet <- as.data.frame(fread("Res/trainTestSplit/training_set.csv"))
testSet <- as.data.frame(fread("Res/trainTestSplit/testing_set.csv"))

trainGroup <- trainSet$group
testGroup <- testSet$group

rownames(trainSet) <- trainSet$V1
rownames(testSet) <- testSet$V1

trainData <- trainSet[,-c(1,2)]
testData <- testSet[,-c(1,2)]

lassoGenes <- fread("Res/LASSO/lassoGenes.txt", header = FALSE)$V1
lassoGenes <- sort(lassoGenes)
biomarkers <- lassoGenes[-c(3,10)]

trainData <- trainData[,which(colnames(trainData) %in% biomarkers)]
testData <- testData[,which(colnames(testData) %in% biomarkers)]

trainData <- trainData[,biomarkers]
testData <- testData[,biomarkers]

trainGroup <- as.factor(trainGroup)
testGroup <- as.factor(testGroup)

levels(trainGroup) <- c("Normal", "Tumor")
levels(testGroup) <- c("Normal", "Tumor")

train <- cbind(group=trainGroup, trainData)
test <- cbind(group=testGroup, testData)
#####################################

#### ROC curve ####
###################
roc.curve <- function(rocs, cols, main, mods) {
  
  png(paste0("Res/trainEvaluateML/roc", main ,".png"), width = 1800, height = 1800, res = 300)
  
  plot.roc(rocs[[1]],
           col=cols[1],
           legacy.axes = FALSE,
           main = main,
           cex.main = 1.6,
           lwd = 2.5,
           cex.lab = 1.3
  )
  
  lines(rocs[[2]], col=cols[2])
  lines(rocs[[3]], col=cols[3])
  
  legend("bottomright",
         legend=c(paste(mods[1], "(AUC:", signif(rocs[[1]]$auc*100, digits = 4), "%)"),
                  paste(mods[2], "(AUC:", signif(rocs[[2]]$auc*100, digits = 4), "%)"),
                  paste(mods[3], "(AUC:", signif(rocs[[3]]$auc*100, digits = 4), "%)")
                  ),
         col=c(cols[1], cols[2], cols[3]),
         lty = 1,
         lwd = 2
  )
  dev.off()
}
###################


#### Random forest ####
#######################

# Tune the random forest
best.mtry <- tuneRF(train[,-1], train[,1], 
                    ntreeTry = 500,
                    stepFactor = 2, 
                    improve = 0.05, 
                    trace = TRUE,
                    plot = FALSE)

best.m <- best.mtry[best.mtry[,2] == min(best.mtry[,2]),1]

best.rf <- randomForest(formula = group ~ ., 
                        data = train, 
                        mtry=best.m[1],
                        importance=TRUE,
                        ntree=500)

# Macke prediction on train data
rf.train.group.true <- train$group
rf.train.group.pred <- as.factor(predict(best.rf, train[,-1]))
rf.train.group.prob <- as.data.frame(predict(best.rf, train[,-1], type = "prob"))

rf.train.roc <- roc(rf.train.group.true ~ as.numeric(rf.train.group.prob$Tumor))

rf.train.cm <- confusionMatrix(rf.train.group.true, rf.train.group.pred, positive = "Tumor")
write.csv(rf.train.cm$table, "Res/trainEvaluateML/rf.train.confusion.matrix.csv")

rf.train.group.true <- as.numeric(rf.train.group.true)
rf.train.group.pred <- as.numeric(rf.train.group.pred)

rf.train.matrics <- data.frame("AUC"=rf.train.roc$auc[1],
                         "accuracy"=Accuracy(rf.train.group.pred, rf.train.group.true),
                         "sensitivity"=Sensitivity(rf.train.group.true, rf.train.group.pred, positive = 2),
                         "specificity"=Specificity(rf.train.group.true, rf.train.group.pred, positive = 2),
                         "precision"=Precision(rf.train.group.true, rf.train.group.pred, positive = 2),
                         "recall"=Recall(rf.train.group.true, rf.train.group.pred, positive = 2),
                         "f1score"=F1_Score(rf.train.group.true, rf.train.group.pred, positive = 2)
                         )

write.csv(rf.train.matrics, "Res/trainEvaluateML/rf.train.matrics.csv")


# Macke prediction on test data
rf.test.group.true <- test$group
rf.test.group.pred <- as.factor(predict(best.rf, test[,-1]))
rf.test.group.prob <- as.data.frame(predict(best.rf, test[,-1], type = "prob"))

rf.test.roc <- roc(rf.test.group.true ~ as.numeric(rf.test.group.prob$Tumor))

rf.test.cm <- confusionMatrix(rf.test.group.true, rf.test.group.pred, positive = "Tumor")
write.csv(rf.test.cm$table, "Res/trainEvaluateML/rf.test.confusion.matrix.csv")

rf.test.group.true <- as.numeric(rf.test.group.true)
rf.test.group.pred <- as.numeric(rf.test.group.pred)

rf.test.matrics <- data.frame("AUC"=rf.test.roc$auc[1],
                              "accuracy"=Accuracy(rf.test.group.pred, rf.test.group.true),
                              "sensitivity"=Sensitivity(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "specificity"=Specificity(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "precision"=Precision(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "recall"=Recall(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "f1score"=F1_Score(rf.test.group.true, rf.test.group.pred, positive = 2)
)

write.csv(rf.test.matrics, "Res/trainEvaluateML/rf.test.matrics.csv")
#######################


#### Support Vector Machine ####
################################

# find optimal cost of misclassification
tune.svm <- tune(svm, group ~ ., 
                 data = train, 
                 kernel = "radial",
                 ranges = list(cost = c(0.001, 0.01, 0.1,1,10,100,1000),
                               gamma = c(0.5,1,2,3,4)),
                 
                 )

svm.fit <- svm(group ~ .,
               data = train,
               kernel = "radial",
               cost = tune.svm$best.parameters$cost,
               gamma = tune.svm$best.parameters$gamma,
               probability = TRUE
               )

# Macke prediction on train data
svm.train.group.true <- train$group
svm.train.group.pred <- unname(predict(svm.fit, newdata = train[,-1]))
svm.train.group.prob <- as.data.frame(attr(predict(svm.fit, newdata = train[,-1], probability = TRUE), "probabilities"))
  
svm.train.roc <- roc(svm.train.group.true ~ as.numeric(svm.train.group.prob$Tumor))

svm.train.cm <- table(svm.train.group.true, svm.train.group.pred)
write.csv(svm.train.cm, "Res/trainEvaluateML/svm.train.confusion.matrix.csv")


svm.train.matrics <- data.frame("AUC"=svm.train.roc$auc[1],
                               "accuracy"=Accuracy(svm.train.group.pred, svm.train.group.true),
                               "sensitivity"=Sensitivity(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "precision"=Precision(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "recall"=Recall(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(svm.train.group.true, svm.train.group.pred, positive = "Tumor")
)

write.csv(svm.train.matrics, "Res/trainEvaluateML/svm.train.matrics.csv")


# Macke prediction on test data
svm.test.group.true <- test$group
svm.test.group.pred <- unname(predict(tune.svm$best.model, newdata = test[,-1]))
svm.test.group.prob <- as.data.frame(attr(predict(svm.fit, newdata = test[,-1], probability = TRUE), "probabilities"))

svm.test.roc <- roc(svm.test.group.true ~ as.numeric(svm.test.group.prob$Tumor))

svm.test.cm <- table(svm.test.group.true, svm.test.group.pred)
write.csv(svm.test.cm, "Res/trainEvaluateML/svm.test.confusion.matrix.csv")

svm.test.matrics <- data.frame("AUC"=svm.test.roc$auc[1],
                               "accuracy"=Accuracy(svm.test.group.pred, svm.test.group.true),
                               "sensitivity"=Sensitivity(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "precision"=Precision(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "recall"=Recall(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(svm.test.group.true, svm.test.group.pred, positive = "Tumor")
)

write.csv(svm.test.matrics, "Res/trainEvaluateML/svm.test.matrics.csv")
################################


#### ROC plot ####
##################
#cols <- c("red", "blue")
#mods <- c("RF", "SVM")

#rocs.train <- list(rf.train.roc, svm.train.roc)
#roc.curve(rocs.train, cols, "Train group", mods)

#rocs.test <- list(rf.test.roc, svm.test.roc)
#roc.curve(rocs.test, cols, "Test group", mods)

saveRDS(best.rf, "Res/trainEvaluateML/RF.rds")
saveRDS(svm.fit, "Res/trainEvaluateML/SVM.rds")
##################