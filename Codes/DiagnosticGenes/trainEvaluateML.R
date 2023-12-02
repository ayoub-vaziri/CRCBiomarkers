library(data.table)
library(ggplot2)
library(caret)
library(pROC)
library(MLmetrics)
library(randomForest)
library(e1071)
library(dplyr)

# Set the current working directory to the project path
setwd(project_path)

set.seed(123)

#### Load and process train data ####
#####################################
trainSet <- as.data.frame(fread("Results/DataProcessing/trainvalidSplit/training_set.csv"))
validSet <- as.data.frame(fread("Results/DataProcessing/trainvalidSplit/validation_set.csv"))

trainGroup <- trainSet$group
validGroup <- validSet$group

rownames(trainSet) <- trainSet$V1
rownames(validSet) <- validSet$V1

trainData <- trainSet[,-c(1,2)]
validData <- validSet[,-c(1,2)]

lassoGenes <- fread("Results/DiagnosticGenes/LASSO/lassoGenes.txt", header = FALSE)$V1
lassoGenes <- sort(lassoGenes)
biomarkers <- lassoGenes[-c(3,10)]

trainData <- trainData[,which(colnames(trainData) %in% biomarkers)]
validData <- validData[,which(colnames(validData) %in% biomarkers)]

trainData <- trainData[,biomarkers]
validData <- validData[,biomarkers]

trainGroup <- as.factor(trainGroup)
validGroup <- as.factor(validGroup)

levels(trainGroup) <- c("Normal", "Tumor")
levels(validGroup) <- c("Normal", "Tumor")

train <- cbind(group=trainGroup, trainData)
valid <- cbind(group=validGroup, validData)
#####################################

#### ROC curve ####
###################
roc.curve <- function(rocs, cols, main, mods) {
  
  png(paste0("Results/DiagnosticGenes/trainEvaluateML/roc", main ,".png"), width = 1800, height = 1800, res = 300)
  
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
write.csv(rf.train.cm$table, "Results/DiagnosticGenes/trainEvaluateML/rf.train.confusion.matrix.csv")

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

write.csv(rf.train.matrics, "Results/DiagnosticGenes/trainEvaluateML/rf.train.matrics.csv")


# Macke prediction on valid data
rf.valid.group.true <- valid$group
rf.valid.group.pred <- as.factor(predict(best.rf, valid[,-1]))
rf.valid.group.prob <- as.data.frame(predict(best.rf, valid[,-1], type = "prob"))

rf.valid.roc <- roc(rf.valid.group.true ~ as.numeric(rf.valid.group.prob$Tumor))

rf.valid.cm <- confusionMatrix(rf.valid.group.true, rf.valid.group.pred, positive = "Tumor")
write.csv(rf.valid.cm$table, "Results/DiagnosticGenes/trainEvaluateML/rf.valid.confusion.matrix.csv")

rf.valid.group.true <- as.numeric(rf.valid.group.true)
rf.valid.group.pred <- as.numeric(rf.valid.group.pred)

rf.valid.matrics <- data.frame("AUC"=rf.valid.roc$auc[1],
                              "accuracy"=Accuracy(rf.valid.group.pred, rf.valid.group.true),
                              "sensitivity"=Sensitivity(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                              "specificity"=Specificity(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                              "precision"=Precision(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                              "recall"=Recall(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                              "f1score"=F1_Score(rf.valid.group.true, rf.valid.group.pred, positive = 2)
)

write.csv(rf.valid.matrics, "Results/DiagnosticGenes/trainEvaluateML/rf.valid.matrics.csv")
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
write.csv(svm.train.cm, "Results/DiagnosticGenes/trainEvaluateML/svm.train.confusion.matrix.csv")


svm.train.matrics <- data.frame("AUC"=svm.train.roc$auc[1],
                               "accuracy"=Accuracy(svm.train.group.pred, svm.train.group.true),
                               "sensitivity"=Sensitivity(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "precision"=Precision(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "recall"=Recall(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(svm.train.group.true, svm.train.group.pred, positive = "Tumor")
)

write.csv(svm.train.matrics, "Results/DiagnosticGenes/trainEvaluateML/svm.train.matrics.csv")


# Macke prediction on valid data
svm.valid.group.true <- valid$group
svm.valid.group.pred <- unname(predict(tune.svm$best.model, newdata = valid[,-1]))
svm.valid.group.prob <- as.data.frame(attr(predict(svm.fit, newdata = valid[,-1], probability = TRUE), "probabilities"))

svm.valid.roc <- roc(svm.valid.group.true ~ as.numeric(svm.valid.group.prob$Tumor))

svm.valid.cm <- table(svm.valid.group.true, svm.valid.group.pred)
write.csv(svm.valid.cm, "Results/DiagnosticGenes/trainEvaluateML/svm.valid.confusion.matrix.csv")

svm.valid.matrics <- data.frame("AUC"=svm.valid.roc$auc[1],
                               "accuracy"=Accuracy(svm.valid.group.pred, svm.valid.group.true),
                               "sensitivity"=Sensitivity(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                               "precision"=Precision(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                               "recall"=Recall(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor")
)

write.csv(svm.valid.matrics, "Results/DiagnosticGenes/trainEvaluateML/svm.valid.matrics.csv")
################################


#### Save trained models ####
#############################
saveRDS(best.rf, "Results/DiagnosticGenes/trainEvaluateML/RF.rds")
saveRDS(svm.fit, "Results/DiagnosticGenes/trainEvaluateML/SVM.rds")
#############################