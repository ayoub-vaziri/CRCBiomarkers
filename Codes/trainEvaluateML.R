library(data.table)
library(ggplot2)
library(caret)
library(pROC)
library(MLmetrics)
library(randomForest)
library(e1071)
library(neuralnet)
library(nnet)
library(dplyr)
library(NeuralNetTools)
library(h2o)

# Set the current working directory to the project path
setwd(project_path)

set.seed(123)

#### Load and process test dataset ####
#######################################
gse21815 <- getGEO("GSE21815", destdir = "Data/")
gse106582 <- getGEO("GSE106582", destdir = "Data/")

gse21815 <- gse21815[[1]]
gse106582 <- gse106582[[1]]

expr21815 <- exprs(gse21815)
expr106582 <- exprs(gse106582)

feat21815 <- fData(gse21815)
feat106582 <- fData(gse106582)

phen21815 <- pData(gse21815)
phen106582 <- pData(gse106582)

expr21815 <- cbind(ID=rownames(expr21815), expr21815)
expr106582 <- cbind(ID=rownames(expr106582), expr106582)

feat21815 <- cbind(ID=feat21815$ID, symbol=feat21815$GENE_SYMBOL)
feat106582 <- cbind(ID=feat106582$ID, symbol=feat106582$Symbol)

expr21815 <- merge(feat21815, expr21815, by="ID")
expr106582 <- merge(feat106582, expr106582, by="ID")

rownames(expr21815) <- expr21815$ID
rownames(expr106582) <- expr106582$ID

expr21815 <- expr21815[,-1]
expr106582 <- expr106582[,-1]

lassoGenes <- fread("Res/LASSO/lassoGenes.txt", header = FALSE)$V1
lassoGenes <- sort(lassoGenes)
biomarkers <- lassoGenes[-c(3,5)]

expr21815 <- expr21815[which(expr21815$symbol %in% biomarkers),]
expr106582 <- expr106582[which(expr106582$symbol %in% biomarkers),]

rownames(expr21815) <- NULL
rownames(expr106582) <- NULL

expr21815 <- expr21815[!duplicated(expr21815$symbol),]
expr106582 <- expr106582[!duplicated(expr106582$symbol),]

rownames(expr21815) <- expr21815$symbol
rownames(expr106582) <- expr106582$symbol

expr21815 <- expr21815[,-1]
expr106582 <- expr106582[,-1]

expr21815 <- as.data.frame(t(expr21815))
expr106582 <- as.data.frame(t(expr106582))

expr21815 <- expr21815[,biomarkers]
expr106582 <- expr106582[,biomarkers]

for (j in 1:ncol(expr21815)) expr21815[,j] <- as.numeric(expr21815[,j])
for (j in 1:ncol(expr106582)) expr106582[,j] <- as.numeric(expr106582[,j])

expr21815 <- na.omit(expr21815)
expr106582 <- na.omit(expr106582)

expr21815 <- cbind(sample=rownames(expr21815), expr21815)
expr106582 <- cbind(sample=rownames(expr106582), expr106582)

phen21815 <- cbind(sample=rownames(phen21815), phen21815)
phen106582 <- cbind(sample=rownames(phen106582), phen106582)

phen21815 <- cbind(sample=phen21815[,"sample"], group=phen21815[,9]) %>% as.data.frame()
phen106582 <- cbind(sample=phen106582[,"sample"], group=phen106582[,37]) %>% as.data.frame()

phen21815$group <- gsub("colorectal cancer", "Tumor", phen21815$group)
phen21815$group <- gsub("normal epithelium", "Normal", phen21815$group)

phen106582$group <- gsub("tumor", "Tumor", phen106582$group)
phen106582$group <- gsub("mucosa", "Normal", phen106582$group)

expr21815 <- merge(phen21815, expr21815, by="sample")
expr106582 <- merge(phen106582, expr106582, by="sample")

rownames(expr21815) <- expr21815$sample
rownames(expr106582) <- expr106582$sample

expr21815 <- expr21815[,-1]
expr21815[,-1] <- log2(expr21815[,-1]+1)
expr106582 <- expr106582[,-1]

write.csv(expr21815, "Results/trainEvaluateML/GSE21815.csv")
write.csv(expr106582, "Results/trainEvaluateML/GSE106582.csv")
#######################################

#### Load and process train and test data ####
##############################################
trainSet <- as.data.frame(fread("Results/batchCorrection/corrected_training_set.csv"))
validSet <- as.data.frame(fread("Results/batchCorrection/corrected_validation_set.csv"))
gse21815 <- as.data.frame(fread("Results/trainEvaluateML/GSE21815.csv"))
gse106582 <- as.data.frame(fread("Results/trainEvaluateML/GSE106582.csv"))

trainGroup <- trainSet$group
validGroup <- validSet$group
gse21815Group <- gse21815$group
gse106582Group <- gse106582$group

rownames(trainSet) <- trainSet$sample
rownames(validSet) <- validSet$sample
rownames(gse21815) <- gse21815$V1
rownames(gse106582) <- gse106582$V1

trainData <- trainSet[,-c(1:4)]
validData <- validSet[,-c(1:4)]
gse21815Data <- gse21815[,-c(1,2)]
gse106582Data <- gse106582[,-c(1,2)]

biomarkers <- fread("Results/LASSO/lassoGenes.txt", header = FALSE)$V1
biomarkers <- sort(biomarkers)
biomarkers <- biomarkers[-c(3,5)]

trainData <- trainData[,which(colnames(trainData) %in% biomarkers)]
validData <- validData[,which(colnames(validData) %in% biomarkers)]
gse21815Data <- gse21815Data[,which(colnames(gse21815Data) %in% biomarkers)]
gse106582Data <- gse106582Data[,which(colnames(gse106582Data) %in% biomarkers)]

trainData <- trainData[,biomarkers]
validData <- validData[,biomarkers]
gse21815Data <- gse21815Data[,biomarkers]
gse106582Data <- gse106582Data[,biomarkers]

trainGroup <- as.factor(trainGroup)
validGroup <- as.factor(validGroup)
gse21815Group <- as.factor(gse21815Group)
gse106582Group <- as.factor(gse106582Group)

levels(trainGroup) <- c("Normal", "Tumor")
levels(validGroup) <- c("Normal", "Tumor")
levels(gse21815Group) <- c("Normal", "Tumor")
levels(gse106582Group) <- c("Normal", "Tumor")

train <- cbind(group=trainGroup, trainData)
valid <- cbind(group=validGroup, validData)
gse21815 <- cbind(group=gse21815Group, gse21815Data)
gse106582 <- cbind(group=gse106582Group, gse106582Data)

test <- rbind(gse106582)
accession <- "GSE106582"

write.csv(train, "Results/trainEvaluateML/train_data.csv")
write.csv(valid, "Results/trainEvaluateML/valid_data.csv")
write.csv(test, "Results/trainEvaluateML/test_data.csv")
##############################################

#### ROC curve ####
###################
roc.curve <- function(rocs, cols, main, mods) {
  
  png(paste0("Results/trainEvaluateML/roc", main ,".png"), width = 1800, height = 1800, res = 300)
  
  plot.roc(rocs[[1]],
           col=cols[1],
           legacy.axes = FALSE,
           main = main,
           cex.main = 1.5,
           lwd = 2.5,
           cex.lab = 1.5
  )
  
  lines(rocs[[2]], col=cols[2])
  lines(rocs[[3]], col=cols[3])
  lines(rocs[[4]], col=cols[4])
  
  legend("bottomright",
         legend=c(paste0(mods[1], " (AUC=", signif(rocs[[1]]$auc*100, digits = 4), "%)"),
                  paste0(mods[2], " (AUC=", signif(rocs[[2]]$auc*100, digits = 4), "%)"),
                  paste0(mods[3], " (AUC=", signif(rocs[[3]]$auc*100, digits = 4), "%)"),
                  paste0(mods[4], " (AUC=", signif(rocs[[4]]$auc*100, digits = 4), "%)")
                  ),
         col=c(cols[1], cols[2], cols[3], cols[4]),
         lty = 1,
         lwd = 3,
         cex = 1.2
         # box.lwd = 0.1,
         # box.lty = 2
  )
  dev.off()
}
###################


#### Random forest ####
#######################
ntree = 500

tuned_rf <- tuneRF(train[,-1], 
                   train$group,
                   ntreeTry = ntree,
                   stepFactor = 1.5,
                   improve = 0.05,
                   trace = TRUE,
                   plot = FALSE
                   )

best.m <- tuned_rf[tuned_rf[,2] == min(tuned_rf[,2]),1][1]

best.rf <- randomForest(formula = group ~ ., 
                        data = train, 
                        mtry=best.m,
                        importance=TRUE,
                        ntree=ntree
                        )

# Macke prediction on train data
rf.train.group.true <- train$group
rf.train.group.pred <- as.factor(predict(best.rf, train[,-1]))
rf.train.group.prob <- as.data.frame(predict(best.rf, train[,-1], type = "prob"))

rf.train.roc <- roc(rf.train.group.true ~ as.numeric(rf.train.group.prob$Tumor), ci = TRUE)

rf.train.cm <- confusionMatrix(rf.train.group.true, rf.train.group.pred, positive = "Tumor")
write.csv(rf.train.cm$table, "Results/trainEvaluateML/rf.train.confusion.matrix.csv")

rf.train.group.true <- as.numeric(rf.train.group.true)
rf.train.group.pred <- as.numeric(rf.train.group.pred)

rf.train.matrics <- data.frame("AUC"=rf.train.roc$auc[1],
                               "CI"=paste0(round(rf.train.roc$ci[1], 4), "-", round(rf.train.roc$ci[3], 4)),
                         "accuracy"=Accuracy(rf.train.group.pred, rf.train.group.true),
                         "sensitivity"=Sensitivity(rf.train.group.true, rf.train.group.pred, positive = 2),
                         "specificity"=Specificity(rf.train.group.true, rf.train.group.pred, positive = 2),
                         "precision"=Precision(rf.train.group.true, rf.train.group.pred, positive = 2),
                         "recall"=Recall(rf.train.group.true, rf.train.group.pred, positive = 2),
                         "f1score"=F1_Score(rf.train.group.true, rf.train.group.pred, positive = 2)
                         )

write.csv(rf.train.matrics, "Results/trainEvaluateML/rf.train.matrics.csv")


# Macke prediction on valid data
rf.valid.group.true <- valid$group
rf.valid.group.pred <- as.factor(predict(best.rf, valid[,-1]))
rf.valid.group.prob <- as.data.frame(predict(best.rf, valid[,-1], type = "prob"))

rf.valid.roc <- roc(rf.valid.group.true ~ as.numeric(rf.valid.group.prob$Tumor), ci = TRUE)

rf.valid.cm <- confusionMatrix(rf.valid.group.true, rf.valid.group.pred, positive = "Tumor")
write.csv(rf.valid.cm$table, "Results/trainEvaluateML/rf.valid.confusion.matrix.csv")

rf.valid.group.true <- as.numeric(rf.valid.group.true)
rf.valid.group.pred <- as.numeric(rf.valid.group.pred)

rf.valid.matrics <- data.frame("AUC"=rf.valid.roc$auc[1],
                               "CI"=paste0(round(rf.valid.roc$ci[1], 4), "-", round(rf.valid.roc$ci[3], 4)),
                              "accuracy"=Accuracy(rf.valid.group.pred, rf.valid.group.true),
                              "sensitivity"=Sensitivity(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                              "specificity"=Specificity(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                              "precision"=Precision(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                              "recall"=Recall(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                              "f1score"=F1_Score(rf.valid.group.true, rf.valid.group.pred, positive = 2)
)

write.csv(rf.valid.matrics, "Results/trainEvaluateML/rf.valid.matrics.csv")

# Macke prediction on test data
rf.test.group.true <- test$group
rf.test.group.pred <- as.factor(predict(best.rf, test[,-1]))
rf.test.group.prob <- as.data.frame(predict(best.rf, test[,-1], type = "prob"))

rf.test.roc <- roc(rf.test.group.true ~ as.numeric(rf.test.group.prob$Tumor), ci = TRUE)

rf.test.cm <- confusionMatrix(rf.test.group.true, rf.test.group.pred, positive = "Tumor")
write.csv(rf.test.cm$table, paste0("Results/trainEvaluateML/rf.test.confusion.matrix.", accession, ".csv"))

rf.test.group.true <- as.numeric(rf.test.group.true)
rf.test.group.pred <- as.numeric(rf.test.group.pred)

rf.test.matrics <- data.frame("AUC"=rf.test.roc$auc[1],
                              "CI"=paste0(round(rf.test.roc$ci[1], 4), "-", round(rf.test.roc$ci[3], 4)),
                              "accuracy"=Accuracy(rf.test.group.pred, rf.test.group.true),
                              "sensitivity"=Sensitivity(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "specificity"=Specificity(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "precision"=Precision(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "recall"=Recall(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "f1score"=F1_Score(rf.test.group.true, rf.test.group.pred, positive = 2)
)

write.csv(rf.test.matrics, paste0("Results/trainEvaluateML/rf.test.matrics.", accession, ".csv"))
#######################


#### Support Vector Machine ####
################################
tune.svm <- tune(svm, group ~ ., 
                 data = train, 
                 kernel = "linear",
                 ranges = list(cost = c(0.1, 1, 10),
                               gamma = c(0.1, 0.01, 0.001))
                 )

svm.fit <- svm(group ~ .,
               data = train,
               kernel = "linear",
               cost = tune.svm$best.parameters$cost,
               gamma = tune.svm$best.parameters$gamma,
               probability = TRUE
               )

# Macke prediction on train data
svm.train.group.true <- train$group
svm.train.group.pred <- unname(predict(svm.fit, newdata = train[,-1]))
svm.train.group.prob <- as.data.frame(attr(predict(svm.fit, newdata = train[,-1], probability = TRUE), "probabilities"))
  
svm.train.roc <- roc(svm.train.group.true ~ as.numeric(svm.train.group.prob$Tumor), ci = TRUE)

svm.train.cm <- table(svm.train.group.true, svm.train.group.pred)
write.csv(svm.train.cm, "Results/trainEvaluateML/svm.train.confusion.matrix.csv")


svm.train.matrics <- data.frame("AUC"=svm.train.roc$auc[1],
                                "CI"=paste0(round(svm.train.roc$ci[1], 4), "-", round(svm.train.roc$ci[3], 4)),
                               "accuracy"=Accuracy(svm.train.group.pred, svm.train.group.true),
                               "sensitivity"=Sensitivity(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "precision"=Precision(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "recall"=Recall(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(svm.train.group.true, svm.train.group.pred, positive = "Tumor")
)

write.csv(svm.train.matrics, "Results/trainEvaluateML/svm.train.matrics.csv")


# Macke prediction on valid data
svm.valid.group.true <- valid$group
svm.valid.group.pred <- unname(predict(tune.svm$best.model, newdata = valid[,-1]))
svm.valid.group.prob <- as.data.frame(attr(predict(svm.fit, newdata = valid[,-1], probability = TRUE), "probabilities"))

svm.valid.roc <- roc(svm.valid.group.true ~ as.numeric(svm.valid.group.prob$Tumor), ci=TRUE)

svm.valid.cm <- table(svm.valid.group.true, svm.valid.group.pred)
write.csv(svm.valid.cm, "Results/trainEvaluateML/svm.valid.confusion.matrix.csv")

svm.valid.matrics <- data.frame("AUC"=svm.valid.roc$auc[1],
                                "CI"=paste0(round(svm.valid.roc$ci[1], 4), "-", round(svm.valid.roc$ci[3], 4)),
                               "accuracy"=Accuracy(svm.valid.group.pred, svm.valid.group.true),
                               "sensitivity"=Sensitivity(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                               "precision"=Precision(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                               "recall"=Recall(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor")
)

write.csv(svm.valid.matrics, "Results/trainEvaluateML/svm.valid.matrics.csv")

# Macke prediction on test data
svm.test.group.true <- test$group
svm.test.group.pred <- unname(predict(tune.svm$best.model, newdata = test[,-1]))
svm.test.group.prob <- as.data.frame(attr(predict(svm.fit, newdata = test[,-1], probability = TRUE), "probabilities"))

svm.test.roc <- roc(svm.test.group.true ~ as.numeric(svm.test.group.prob$Tumor), ci=TRUE)

svm.test.cm <- table(svm.test.group.true, svm.test.group.pred)
write.csv(svm.test.cm, paste0("Results/trainEvaluateML/svm.test.confusion.matrix.", accession, ".csv"))

svm.test.matrics <- data.frame("AUC"=svm.test.roc$auc[1],
                               "CI"=paste0(round(svm.test.roc$ci[1], 4), "-", round(svm.test.roc$ci[3], 4)),
                               "accuracy"=Accuracy(svm.test.group.pred, svm.test.group.true),
                               "sensitivity"=Sensitivity(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "precision"=Precision(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "recall"=Recall(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(svm.test.group.true, svm.test.group.pred, positive = "Tumor")
)

write.csv(svm.test.matrics, paste0("Results/trainEvaluateML/svm.test.matrics.", accession, ".csv"))
################################


#### Artificial Neural Network ####
###################################

# Encode as a one hot vector multilabel data
train.onehot <- cbind(class.ind(as.factor(train$group)), train[,-1])
valid.onehot <- cbind(class.ind(as.factor(valid$group)), valid[,-1])
test.onehot <- cbind(class.ind(as.factor(test$group)), test[,-1])

# Set up formula
n <- names(train.onehot)
f <- as.formula(paste("Normal + Tumor ~", paste(n[!n %in% c("Normal", "Tumor")], collapse = " + ")))

nn <- neuralnet(f,
                data = train.onehot,
                hidden = c(5),
                act.fct = "tanh",
                linear.output = FALSE)

png("Results/trainEvaluateML/annplot.png", height = 1600, width = 2200, res = 300)
plotnet(nn,
        cex_val = 0.8,
        circle_cex = 4,
        circle_col = "pink",
        bord_col = "pink",
        neg_col = "grey",
        nid = TRUE
        )
dev.off()

# Macke prediction on train data
nn.train.group.prob <- as.data.frame(predict(nn, train.onehot[,3:ncol(train.onehot)]))
colnames(nn.train.group.prob) <- c("Normal", "Tumor")

nn.train.group.pred <- max.col(nn.train.group.prob) |> as.factor()
nn.train.group.true <- max.col(train.onehot[,1:2]) |> as.factor()

levels(nn.train.group.pred) <- c("Normal", "Tumor")
levels(nn.train.group.true) <- c("Normal", "Tumor")

nn.train.cm <- confusionMatrix(nn.train.group.true, nn.train.group.pred, positive = "Tumor")
write.csv(nn.train.cm$table, "Results/trainEvaluateML/nn.train.confusion.matrix.csv")

nn.train.roc <- roc(nn.train.group.true ~ as.numeric(nn.train.group.prob$Tumor), ci = TRUE)

nn.train.matrics <- data.frame("AUC"=nn.train.roc$auc[1],
                               "CI"=paste0(round(nn.train.roc$ci[1], 4), "-", round(nn.train.roc$ci[3], 4)),
                               "accuracy"=Accuracy(nn.train.group.pred, nn.train.group.true),
                               "sensitivity"=Sensitivity(nn.train.group.true, nn.train.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(nn.train.group.true, nn.train.group.pred, positive = "Tumor"),
                               "precision"=Precision(nn.train.group.true, nn.train.group.pred, positive = "Tumor"),
                               "recall"=Recall(nn.train.group.true, nn.train.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(nn.train.group.true, nn.train.group.pred, positive = "Tumor")
)

write.csv(nn.train.matrics, "Results/trainEvaluateML/nn.train.matrics.csv")


# Macke prediction on valid data
nn.valid.group.prob <- as.data.frame(predict(nn, valid.onehot[,3:ncol(valid.onehot)]))
colnames(nn.valid.group.prob) <- c("Normal", "Tumor")

nn.valid.group.pred <- max.col(nn.valid.group.prob) |> as.factor()
nn.valid.group.true <- max.col(valid.onehot[,1:2]) |> as.factor()

levels(nn.valid.group.pred) <- c("Normal", "Tumor")
levels(nn.valid.group.true) <- c("Normal", "Tumor")

nn.valid.cm <- confusionMatrix(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor")
write.csv(nn.valid.cm$table, "Results/trainEvaluateML/nn.valid.confusion.matrix.csv")

nn.valid.roc <- roc(nn.valid.group.true ~ as.numeric(nn.valid.group.prob$Tumor), ci = TRUE)

nn.valid.matrics <- data.frame("AUC"=nn.valid.roc$auc[1],
                               "CI"=paste0(round(nn.valid.roc$ci[1], 4), "-", round(nn.valid.roc$ci[3], 4)),
                              "accuracy"=Accuracy(nn.valid.group.pred, nn.valid.group.true),
                              "sensitivity"=Sensitivity(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor"),
                              "specificity"=Specificity(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor"),
                              "precision"=Precision(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor"),
                              "recall"=Recall(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor"),
                              "f1score"=F1_Score(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor")
)

write.csv(nn.valid.matrics, "Results/trainEvaluateML/nn.valid.matrics.csv")

# Macke prediction on test data
nn.test.group.prob <- as.data.frame(predict(nn, test.onehot[,3:ncol(test.onehot)]))
colnames(nn.test.group.prob) <- c("Normal", "Tumor")

nn.test.group.pred <- max.col(nn.test.group.prob) |> as.factor()
nn.test.group.true <- max.col(test.onehot[,1:2]) |> as.factor()

levels(nn.test.group.pred) <- c("Normal", "Tumor")
levels(nn.test.group.true) <- c("Normal", "Tumor")

nn.test.cm <- confusionMatrix(nn.test.group.true, nn.test.group.pred, positive = "Tumor")
write.csv(nn.test.cm$table, paste0("Results/trainEvaluateML/nn.test.confusion.matrix.", accession, ".csv"))

nn.test.roc <- roc(nn.test.group.true ~ as.numeric(nn.test.group.prob$Tumor), ci = TRUE)

nn.test.matrics <- data.frame("AUC"=nn.test.roc$auc[1],
                              "CI"=paste0(round(nn.test.roc$ci[1], 4), "-", round(nn.test.roc$ci[3], 4)),
                              "accuracy"=Accuracy(nn.test.group.pred, nn.test.group.true),
                              "sensitivity"=Sensitivity(nn.test.group.true, nn.test.group.pred, positive = "Tumor"),
                              "specificity"=Specificity(nn.test.group.true, nn.test.group.pred, positive = "Tumor"),
                              "precision"=Precision(nn.test.group.true, nn.test.group.pred, positive = "Tumor"),
                              "recall"=Recall(nn.test.group.true, nn.test.group.pred, positive = "Tumor"),
                              "f1score"=F1_Score(nn.test.group.true, nn.test.group.pred, positive = "Tumor")
)

write.csv(nn.test.matrics, paste0("Results/trainEvaluateML/nn.test.matrics.", accession, ".csv"))
###################################

#### Gradient Boosting Machine ####
###################################
h2o.init()

train <- h2o.importFile("Results/trainEvaluateML/train_data.csv")
valid <- h2o.importFile("Results/trainEvaluateML/valid_data.csv")
test <- h2o.importFile("Results/trainEvaluateML/test_data.csv")

response <- "group"
predictors <- colnames(train)[!colnames(train) %in% c("C1", response)]

train[,response] <- as.factor(train[,response])
valid[,response] <- as.factor(valid[,response])
test[,response] <- as.factor(test[,response])

gbm <- h2o.gbm(x = predictors,
               y = response,
               training_frame = train,
               ntrees = 100
)

# Macke prediction on train data
gbm.train.group.true <- as.data.frame(train$group)$group
gbm.train.group.pred <- as.data.frame(h2o.predict(object = gbm, newdata = train[,-c(1,2)]))$predict
gbm.train.group.prob <- as.data.frame(h2o.predict(object = gbm, newdata = train[,-c(1,2)]))

gbm.train.roc <- roc(gbm.train.group.true ~ as.numeric(gbm.train.group.prob$Tumor), ci = TRUE)

gbm.train.cm <- confusionMatrix(gbm.train.group.true, gbm.train.group.pred, positive = "Tumor")
write.csv(gbm.train.cm$table, "Results/trainEvaluateML/gbm.train.confusion.matrix.csv")

gbm.train.group.true <- as.numeric(gbm.train.group.true)
gbm.train.group.pred <- as.numeric(gbm.train.group.pred)

gbm.train.matrics <- data.frame("AUC"=gbm.train.roc$auc[1],
                                "CI"=paste0(round(gbm.train.roc$ci[1], 4), "-", round(gbm.train.roc$ci[3], 4)),
                                "accuracy"=Accuracy(gbm.train.group.pred, gbm.train.group.true),
                                "sensitivity"=Sensitivity(gbm.train.group.true, gbm.train.group.pred, positive = 2),
                                "specificity"=Specificity(gbm.train.group.true, gbm.train.group.pred, positive = 2),
                                "precision"=Precision(gbm.train.group.true, gbm.train.group.pred, positive = 2),
                                "recall"=Recall(gbm.train.group.true, gbm.train.group.pred, positive = 2),
                                "f1score"=F1_Score(gbm.train.group.true, gbm.train.group.pred, positive = 2)
)

write.csv(gbm.train.matrics, "Results/trainEvaluateML/gbm.train.matrics.csv")

# Macke prediction on valid data
gbm.valid.group.true <- as.data.frame(valid$group)$group
gbm.valid.group.pred <- as.data.frame(h2o.predict(object = gbm, newdata = valid[,-c(1,2)]))$predict
gbm.valid.group.prob <- as.data.frame(h2o.predict(object = gbm, newdata = valid[,-c(1,2)]))

gbm.valid.roc <- roc(gbm.valid.group.true ~ as.numeric(gbm.valid.group.prob$Tumor), ci = TRUE)

gbm.valid.cm <- confusionMatrix(gbm.valid.group.true, gbm.valid.group.pred, positive = "Tumor")
write.csv(gbm.valid.cm$table, "Results/trainEvaluateML/gbm.valid.confusion.matrix.csv")

gbm.valid.group.true <- as.numeric(gbm.valid.group.true)
gbm.valid.group.pred <- as.numeric(gbm.valid.group.pred)

gbm.valid.matrics <- data.frame("AUC"=gbm.valid.roc$auc[1],
                                "CI"=paste0(round(gbm.valid.roc$ci[1], 4), "-", round(gbm.valid.roc$ci[3], 4)),
                                "accuracy"=Accuracy(gbm.valid.group.pred, gbm.valid.group.true),
                                "sensitivity"=Sensitivity(gbm.valid.group.true, gbm.valid.group.pred, positive = 2),
                                "specificity"=Specificity(gbm.valid.group.true, gbm.valid.group.pred, positive = 2),
                                "precision"=Precision(gbm.valid.group.true, gbm.valid.group.pred, positive = 2),
                                "recall"=Recall(gbm.valid.group.true, gbm.valid.group.pred, positive = 2),
                                "f1score"=F1_Score(gbm.valid.group.true, gbm.valid.group.pred, positive = 2)
)

write.csv(gbm.valid.matrics, "Results/trainEvaluateML/gbm.valid.matrics.csv")

# Macke prediction on test data
gbm.test.group.true <- as.data.frame(test$group)$group
gbm.test.group.pred <- as.data.frame(h2o.predict(object = gbm, newdata = test[,-c(1,2)]))$predict
gbm.test.group.prob <- as.data.frame(h2o.predict(object = gbm, newdata = test[,-c(1,2)]))

gbm.test.roc <- roc(gbm.test.group.true ~ as.numeric(gbm.test.group.prob$Tumor), ci = TRUE)

gbm.test.cm <- confusionMatrix(gbm.test.group.true, gbm.test.group.pred, positive = "Tumor")
write.csv(gbm.test.cm$table, paste0("Results/trainEvaluateML/gbm.test.confusion.matrix.", accession, ".csv"))

gbm.test.group.true <- as.numeric(gbm.test.group.true)
gbm.test.group.pred <- as.numeric(gbm.test.group.pred)

gbm.test.matrics <- data.frame("AUC"=gbm.test.roc$auc[1],
                               "CI"=paste0(round(gbm.test.roc$ci[1], 4), "-", round(gbm.test.roc$ci[3], 4)),
                               "accuracy"=Accuracy(gbm.test.group.pred, gbm.test.group.true),
                               "sensitivity"=Sensitivity(gbm.test.group.true, gbm.test.group.pred, positive = 2),
                               "specificity"=Specificity(gbm.test.group.true, gbm.test.group.pred, positive = 2),
                               "precision"=Precision(gbm.test.group.true, gbm.test.group.pred, positive = 2),
                               "recall"=Recall(gbm.test.group.true, gbm.test.group.pred, positive = 2),
                               "f1score"=F1_Score(gbm.test.group.true, gbm.test.group.pred, positive = 2)
)

write.csv(gbm.test.matrics, paste0("Results/trainEvaluateML/gbm.test.matrics.", accession, ".csv"))
###################################


#### ROC plot ####
##################
cols <- c("blue", "green", "red", "black")
mods <- c("RF", "SVM", "ANN", "GBM")

rocs.train <- list(rf.train.roc, svm.train.roc, nn.train.roc, gbm.train.roc)
roc.curve(rocs.train, cols, "Training Dataset", mods)

rocs.valid <- list(rf.valid.roc, svm.valid.roc, nn.valid.roc, gbm.valid.roc)
roc.curve(rocs.valid, cols, "Validation Dataset", mods)

rocs.test <- list(rf.test.roc, svm.test.roc, nn.test.roc, gbm.test.roc)
roc.curve(rocs.test, cols, paste0("Test Dataset (", accession, ")"), mods)

saveRDS(best.rf, "Results/trainEvaluateML/RF.rds")
saveRDS(svm.fit, "Results/trainEvaluateML/SVM.rds")
saveRDS(nn, "Results/trainEvaluateML/ANN.rds")
saveRDS(gbm, "Results/trainEvaluateML/GBM.rds")
##################
