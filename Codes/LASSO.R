library(data.table)
library(glmnet)
library(forcats)
library(ggplot2)
library(dplyr)

# Set the current working directory to the project path
setwd(project_path)

#### Load and process train data ####
#####################################
trainSet <- as.data.frame(fread("Results/batchCorrection/corrected_training_set.csv"))

trainGroup <- trainSet$group

rownames(trainSet) <- trainSet$sample

trainData <- trainSet[,-c(1:4)]

# essentialNodes <- fread("Res/PPI/essentialNodes.txt", header = FALSE)$V1

keyGenes <- fread("Results/centralityAnalysis/essentialNodes.txt", header = FALSE)$V1

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

png(filename = "Results/LASSO/Lambda_BinomialDeviance.png",
    width = 1800, height = 1800, res = 300)
plot(cv.lasso, axes=T, cex.lab = 1.4, cex.axis = 1.2)
dev.off()

write.table(cv.lasso$lambda.min, "Results/LASSO/lambda_min.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Final model with lambda.min
lasso.model <- glmnet(x, as.factor(y),
                      alpha = 1,
                      family = "binomial",
                      lambda = cv.lasso$lambda.min)


# Non-zero coefficients
beta <- lasso.model$beta[,1]
coef.nonzero <- beta[which(beta != 0)]
df.coef.nonzero <- data.frame(gene=names(coef.nonzero), coefficient=unname(coef.nonzero))
write.csv(df.coef.nonzero, "Results/LASSO/coefficients.csv")

beta.nonzero <- names(which(beta != 0))

write.table(x = beta.nonzero,
            file = "Results/LASSO/lassoGenes.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

gene.coef <- beta[which(beta != 0)]
dat <- data.frame(Gene=names(gene.coef), Beta=unname(gene.coef))
dat <- dat[order(dat$Beta),]
rownames(dat) <- NULL

beta <- as.numeric(union(substr(dat$Beta[1], 1, 6), substr(dat$Beta[c(1:11)], 1, 5)))
beta <- beta[-2]

png("Results/LASSO/nonZeroCoefGenes.png", width = 1600, height = 1400, res = 300)
dat %>% 
  arrange(Beta) %>% 
  mutate(Gene = factor(Gene, levels=dat$Gene)) %>% 
  ggplot(aes(fill=ifelse(Beta > 0, 1, "#1F77B4"), y=Gene, x=Beta)) + 
  # fill=ifelse(Beta > 0, "#1F77B4", 1)
  # fill=Gene
  geom_bar(stat="identity", width = 0.65) +
  geom_text(stat='identity', aes(label=beta), hjust=ifelse(beta>0, -0.1, 1.1), position = position_dodge2(width = 0.8), size=3) +
  theme_bw() +
  xlab("Coefficient") +
  ylab("") +
  scale_x_continuous(limits = c(-0.5, 2.2)) +
  theme(axis.title.x = element_text(size = 13),
        axis.text.x = element_text(colour = "black", angle = 0, hjust = 1, size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        legend.position = "none") 
  # geom_text(aes(label=signif(Beta, 3)), vjust = ifelse(dat$Beta > 0, -0.5, 1.2), size = 2)
dev.off()
##########################