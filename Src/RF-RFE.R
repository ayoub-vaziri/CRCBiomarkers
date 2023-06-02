library(data.table)
library(caret)
library(ggplot2)
library(dplyr)

set.seed(123)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load and process train data ####
#####################################
train_data <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
rownames(train_data) <- train_data$V1
train_data <- train_data[,-1]
top20 <- fread("Res/PPI/top20.txt", header = FALSE)$V1
train_data <- train_data[,top20]

status <- fread("Res/exprBatchBioData/bioColorectal.txt", header = FALSE)$V1
status <- as.factor(status)
levels(status) <- c(0, 1)
status <- as.numeric(as.character(status))

train_data <- cbind(STATUS=status, train_data)

write.csv(train_data, "Res/RF/train_data.csv")
#####################################

#### Random forest ####
#######################
x <- train_data[,-1]
y <- as.factor(train_data[,1])

doMC::registerDoMC(cores = 2)

rfPredictor = rfeControl(functions = rfFuncs, method = "cv", number = 10)

rfeSelector <- rfe(x, y, 
                   sizes = 1:20,
                   rfeControl = rfPredictor
                   )
           
print(rfeSelector)

head(rfeSelector$variables, 20)

png(filename = "Res/RF/Accuracyplot.png", width = 1600, height = 1600, res = 300)
ggplot(rfeSelector) +
  annotate("text", x=rfeSelector$bestSubset, y=0.951, label=paste0("N=",rfeSelector$bestSubset), col="red") +
  geom_line(color="forestgreen") +
  geom_point(color="forestgreen") +
  theme_bw(base_size = 10) +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"))

dev.off()

write.table(x = predictors(rfeSelector),
            file = "Res/RF/rfGenes.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

varimp <- head(varImp(rfeSelector), rfeSelector$optsize)
varimp$var <- rownames(varimp)
varimp <- varimp[order(varimp$Overall),]
rownames(varimp) <- NULL

png("Res/RF/variableImportance.png", width = 1600, height = 1400, res = 300)
varimp %>% 
  arrange(Overall) %>% 
  mutate(var = factor(var, levels=varimp$var)) %>% 
  ggplot(aes(y=Overall, x=var)) + 
  geom_bar(stat="identity", fill = "#1F77B4", alpha = 0.8, width = 0.3) +
  theme_classic() +
  xlab("") +
  ylab("Variable importance") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(colour = "black", angle = 0, hjust = 1, size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        legend.position = "none") +
  coord_flip()
dev.off()
#######################