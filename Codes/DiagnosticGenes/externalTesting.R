library(MLmetrics)
library(pROC)
library(randomForest)
library(caret)
library(data.table)
library(GEOquery)
library(dplyr)
library(preprocessCore)

# Set the current working directory to the project path
setwd(project_path)

#### Load and process external testing dataset ####
###################################################
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
biomarkers <- lassoGenes[-c(3,10)]

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
expr106582 <- expr106582[,-1]
###################################################

#### ML models evaluation ####
##############################
# Load trained RF and SVM models
RF <- readRDS("Results/DiagnosticGenes/trainEvaluateML/RF.rds")
SVM <- readRDS("Results/DiagnosticGenes/trainEvaluateML/SVM.rds")

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
  write.csv(cm, paste0("Results/DiagnosticGenes/externalTesting/", accession, ".ConfusionMatrix.csv"))
  mats <- data.frame("AUC"=roc$auc[1],
                     "accuracy"=Accuracy(y.pred, y.true),
                     "sensitivity"=Sensitivity(y.true, y.pred, positive = "Tumor"),
                     "specificity"=Specificity(y.true, y.pred, positive = "Tumor"),
                     "precision"=Precision(y.true, y.pred, positive = "Tumor"),
                     "recall"=Recall(y.true, y.pred, positive = "Tumor"),
                     "f1score"=F1_Score(y.true, y.pred, positive = "Tumor")
                     )

  write.csv(mats, paste0("Results/DiagnosticGenes/externalTesting/", accession, ".Metrics.csv"))
  return(mats)
}


svm.metrics <- function(x, y, accession, model) {
  y.true <- as.factor(y)
  y.pred <- unname(predict(model, newdata = x))
  y.prob <- as.data.frame(attr(predict(model, newdata = x, probability = TRUE), "probabilities"))
  
  roc <- roc(y.true ~ as.numeric(y.prob$Tumor))
  
  cm <- table(y.true, y.pred)
  write.csv(cm, paste0("Results/DiagnosticGenes/externalTesting/", accession, ".ConfusionMatrix.csv"))
  
  mats <- data.frame(
                    "AUC"=roc$auc[1],
                    "accuracy"=Accuracy(y.pred, y.true),
                    "sensitivity"=Sensitivity(y.true, y.pred, positive = "Tumor"),
                    "specificity"=Specificity(y.true, y.pred, positive = "Tumor"),
                    "precision"=Precision(y.true, y.pred, positive = "Tumor"),
                    "recall"=Recall(y.true, y.pred, positive = "Tumor"),
                    "f1score"=F1_Score(y.true, y.pred, positive = "Tumor")
                    )
  write.csv(mats, paste0("Results/DiagnosticGenes/externalTesting/", accession, ".Metrics.csv"))
  return(mats)
}
##############################


#### Evaluation of ML models on external testing set ####
#########################################################
expr21815[,-1] <- log2trans(expr21815[,-1])
expr106582[,-1] <- log2trans(expr106582[,-1])

mergedData <- rbind(expr21815, expr106582)

x <- mergedData[,-1]
y <- mergedData[,1]

rf.metrics(x, y, "External testing set for RF", RF)
svm.metrics(x, y, "External testing set for SVM", SVM)
#########################################################