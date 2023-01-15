library(limma)
library(readr)
library(data.table)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

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
  write_tsv(up.regulated, "Res/limma/topTable.up.tsv")
  
  down.regulated <- subset(tT, adj.P.Val < 0.05 & logFC < -1)
  write_tsv(down.regulated, "Res/limma/topTable.down.tsv")

  up <- unique(up.regulated$Gene.symbol)
  down <- unique(down.regulated$Gene.symbol)
  updown <- union(up,down)
  
  write.table(up, "Res/limma/up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(down, "Res/limma/down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(updown, "Res/limma/updown.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

exprColorectal <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
bioColorectal <- fread("Res/exprBatchBioData/bioColorectal.txt")$V1

rownames(exprColorectal) <- exprColorectal$V1
exprColorectal <- exprColorectal[,-1]
exprColorectal <- as.data.frame(t(exprColorectal))

bioColorectal <- gsub(1, "N", bioColorectal)
bioColorectal <- gsub(2, "C", bioColorectal)
bioColorectal <- as.factor(bioColorectal)
levels(bioColorectal) <- c("Tumor", "Normal")

exprColorectal <- log2trans(exprColorectal)
tT <- degs(exprColorectal, bioColorectal)

tT$`Gene.symbol` <- rownames(tT)
tT <- tT[, c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]
rownames(tT) <- NULL

write_tsv(tT, "Res/limma/topTable.tsv")

updown(tT)
##########################################