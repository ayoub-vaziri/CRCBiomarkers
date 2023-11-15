suppressMessages(library(randomForest, quietly = T))

setwd("~/CRC/")

# log2 transform
log2trans <- function(expr) {
  quan <- as.numeric(quantile(expr, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (quan[5] > 100) || (quan[6]-quan[1] > 50 && qaun[2] > 0)
  if(LogC) 
    expr <- log2(expr+1) 
  return(expr)
}

# Model prediction
pred <- function(model, x) {
  x <- log2trans(x)
  y.pred <- as.character(predict(model, x))
  y.prob <- as.data.frame(predict(model, x, type = "prob"))
  return(list(y.pred, y.prob))
}

biomarkers <- c("CDC25B", "CDK4", "MMP1", "MMP7", "PLAU", "SLC7A5", "TEAD4", "TRIB3")

exprData <- read.table("exprData")$V1
model <- readRDS("TrainedRF.rds")

x <- as.data.frame(matrix(biomarkers, nrow = 1, ncol = 8))
colnames(x) <- x[1,]
x <- as.data.frame(x[-1,])
x[1,] <- exprData
for(i in 1:8) x[,i] <- as.numeric(x[,i])

y.pred <- pred(model, x)[[1]]
y.prob <- pred(model, x)[[2]]

cat("\n=============================\n")
if(y.pred == "Normal") {
  cat("SAMPLE STATUS REPORT:\n")
  cat("\tNORMAL\n")
  cat("PROBABILITY:\n")
  cat("\tNormal: ", as.numeric(y.prob[1]), "\n\tTumor : ", as.numeric(y.prob[2]))
} else {
  cat("SAMPLE STATUS REPORT:\n")
  cat("\tTUMOR\n")
  cat("PROBABILITY:\n")
  cat("\tNormal: ", as.numeric(y.prob[1]), "\n\tTumor : ", as.numeric(y.prob[2])) 
}
cat("\n=============================\n\n")
