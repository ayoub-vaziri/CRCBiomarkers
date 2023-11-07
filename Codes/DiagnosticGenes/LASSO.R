library(data.table)
library(glmnet)
library(forcats)
library(ggplot2)
library(dplyr)

set.seed(123)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load and process train data ####
#####################################
trainSet <- as.data.frame(fread("Res/trainTestSplit/training_set.csv"))

trainGroup <- trainSet$group

rownames(trainSet) <- trainSet$V1

trainData <- trainSet[,-c(1,2)]

keyGenes <- fread("Res/centralityAnalysis/keyGenes.txt", header = FALSE)$V1

trainData <- trainData[,which(colnames(trainData) %in% keyGenes)]

trainGroup <- as.factor(trainGroup)

levels(trainGroup) <- c("Normal", "Tumor")

train <- cbind(group=trainGroup, trainData)
#####################################


#### LASSO logistic regression ####
###################################
set.seed(123)

x <- model.matrix(group ~ ., train)[,-1]
y <- ifelse(train$group == "Tumor", 1, 0)

cv.lasso <- cv.glmnet(x, y, 
                      alpha = 1,
                      family = "binomial",
                      type.measure = "deviance",
                      nfolds = 10)

png(filename = "Res/LASSO/Lambda_BinomialDeviance.png",
    width = 1800, height = 1800, res = 300)
plot(cv.lasso)
dev.off()

write.table(cv.lasso$lambda.min, "Res/LASSO/lambda_min.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Final model with lambda.min
lasso.model <- glmnet(x, as.factor(y),
                      alpha = 1,
                      family = "binomial",
                      lambda = cv.lasso$lambda.min)


# Non-zero coefficients
beta <- lasso.model$beta[,1]
coef.nonzero <- beta[which(beta != 0)]
df.coef.nonzero <- data.frame(gene=names(coef.nonzero), coefficient=unname(coef.nonzero))
write.csv(df.coef.nonzero, "Res/LASSO/coefficients.csv")

beta.nonzero <- names(which(beta != 0))
write.table(x = beta.nonzero,
            file = "Res/LASSO/lassoGenes.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

gene.coef <- beta[which(beta != 0)]
dat <- data.frame(Gene=names(gene.coef), Beta=unname(gene.coef))
dat <- dat[order(dat$Beta),]
rownames(dat) <- NULL

beta <- as.numeric(union(substr(dat$Beta[1], 1, 6), substr(dat$Beta[c(2:10)], 1, 5)))

png("Res/LASSO/nonZeroCoefGenes.png", width = 1600, height = 1400, res = 300)
dat %>% 
  arrange(Beta) %>% 
  mutate(Gene = factor(Gene, levels=dat$Gene)) %>% 
  ggplot(aes(fill=ifelse(Beta > 0, 1, "#1F77B4"), y=Beta, x=Gene)) + 
  # fill=ifelse(Beta > 0, "#1F77B4", 1)
  # fill=Gene
  geom_bar(stat="identity", width = 0.85) +
  geom_text(stat='identity', aes(label=beta), vjust=ifelse(beta>0,-0.4,1.3), position = position_dodge2(width = 0.6), size=3) +
  theme_bw() +
  xlab("") +
  ylab("Coefficient") +
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        legend.position = "none") 
  # geom_text(aes(label=signif(Beta, 3)), vjust = ifelse(dat$Beta > 0, -0.5, 1.2), size = 2)
dev.off()
##########################