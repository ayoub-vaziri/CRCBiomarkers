library(data.table)
library(glmnet)
library(dplyr)
library(forcats)
library(ggplot2)

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
#####################################

#### LASSO regression ####
##########################
x <- model.matrix(STATUS ~ ., data = train_data)[,-1]
y <- train_data$STATUS

cvfit <- cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 10)

coefs <- as.data.frame(coef(cvfit$glmnet.fit, s = cvfit$lambda.min)[,"s1"])
colnames(coefs) <- "BETA"
coefs <- cbind(GENE=rownames(coefs), coefs)
coefs <- coefs[-1,]
rownames(coefs) <- NULL
nonzeroCoefGenes <- coefs$GENE[which(coefs$BETA != 0)]

pdf("Res/LASSO/Lambda_BinomialDeviance.pdf")
plot(cvfit)
dev.off()

png(filename = "Res/LASSO/Lambda_BinomialDeviance.png",
    width = 1800, height = 1800, res = 300)
plot(cvfit)
dev.off()

write.table(x = nonzeroCoefGenes,
            file = "Res/LASSO/lassoGenes.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

gene <- coefs$GENE[which(coefs$BETA != 0)]
beta <- coefs$BETA[which(coefs$BETA != 0)]
dat <- data.frame(Gene=gene, Beta=beta)
dat <- dat[order(dat$Beta),]

pdf("Res/LASSO/nonZeroCoefGenes.pdf")
dat %>% 
  arrange(Beta) %>% 
  mutate(Gene = factor(Gene, levels=dat$Gene)) %>% 
  ggplot(aes(fill=Gene, y=Beta, x=Gene)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("") +
  ylab("Coefficient") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(colour = "black", size = 10),
        legend.position = "none")
dev.off()

png("Res/LASSO/nonZeroCoefGenes.png", width = 1600, height = 1400, res = 300)
dat %>% 
  arrange(Beta) %>% 
  mutate(Gene = factor(Gene, levels=dat$Gene)) %>% 
  ggplot(aes(fill=Gene, y=Beta, x=Gene)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("") +
  ylab("Coefficient") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(colour = "black", size = 10),
        legend.position = "none")
dev.off()
##########################