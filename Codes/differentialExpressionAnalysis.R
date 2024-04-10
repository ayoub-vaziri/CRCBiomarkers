library(limma)
library(readr)
library(data.table)

# Set the current working directory to the project path
setwd(project_path)

#### Differential expression analysis ####
##########################################

# log2 transform
log2trans <- function(expr) {
  quan <- as.numeric(quantile(expr, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (quan[5] > 100) || (quan[6]-quan[1] > 50 && qaun[2] > 0)
  if(LogC) { 
    expr[which(expr <= 0)] <- NaN
    expr <- log2(expr) 
  }
  return(expr)
}

# differential expression analysis
degs <- function(expr, group) {
  design.mat <- model.matrix(~ group + 0)
  colnames(design.mat) <- levels(group)
  lmfit <- lmFit(expr, design.mat)
  cont.mat <- makeContrasts(contrasts = "Tumor-Normal", levels = design.mat)
  cont.fit <- contrasts.fit(lmfit, cont.mat)
  fit <- eBayes(cont.fit, 0.01)
  toptable <- topTable(fit = fit, number = Inf, adjust.method = "fdr", sort.by="B")
  return(toptable)
}

# Find up & down regulated genes
updown <- function(tT) {
  up.regulated <- subset(tT, adj.P.Val < 0.05 & logFC > 1)
  write_tsv(up.regulated, "Results/differentialExpressionAnalysis/topTable.up.tsv")
  
  down.regulated <- subset(tT, adj.P.Val < 0.05 & logFC < -1)
  write_tsv(down.regulated, "Results/differentialExpressionAnalysis/topTable.down.tsv")

  up <- unique(up.regulated$Gene.symbol)
  down <- unique(down.regulated$Gene.symbol)
  updown <- union(up,down)
  
  write.table(up, "Results/differentialExpressionAnalysis/up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(down, "Results/differentialExpressionAnalysis/down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(updown, "Results/differentialExpressionAnalysis/updown.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

trainingSet <- as.data.frame(fread("Results/batchCorrection/corrected_training_set.csv"))

group <- trainingSet$group
group <- gsub(1, "N", group)
group <- gsub(2, "C", group)
group <- as.factor(group)
levels(group) <- c("Tumor", "Normal")

rownames(trainingSet) <- trainingSet$V1
trainingSet <- trainingSet[,-c(1,2)]
trainingSet <- as.data.frame(t(trainingSet))

trainingSet <- log2trans(trainingSet)
tT <- degs(trainingSet, group)

tT$`Gene.symbol` <- rownames(tT)
tT <- tT[, c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]
rownames(tT) <- NULL

write_tsv(tT, "Results/differentialExpressionAnalysis/topTable.tsv")
write_csv(tT, "Results/differentialExpressionAnalysis/topTable.csv")

updown(tT)
##########################################