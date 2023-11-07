library(data.table)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

samples <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
rownames(samples) <- samples$sample
samples <- samples[,-c(1,2,3)]

cancerSamples <- rownames(samples)[which(samples$group == 2)]
normalSamples <- rownames(samples)[which(samples$group == 1)]

trainingCancerSamples <- cancerSamples[sample(1:length(cancerSamples), ceiling(0.7*length(cancerSamples)))]
testingCancerSamples <- cancerSamples[-which(cancerSamples %in% trainingCancerSamples)]

trainingNormalSamples <- normalSamples[sample(1:length(normalSamples), ceiling(0.7*length(normalSamples)))]
testingNormalSamples <- normalSamples[-which(normalSamples %in% trainingNormalSamples)]

trainingSamples <- union(trainingCancerSamples, trainingNormalSamples)
testingSamples <- union(testingCancerSamples, testingNormalSamples)

trainingSet <- samples[which(rownames(samples) %in% trainingSamples),]
testingSet <- samples[which(rownames(samples) %in% testingSamples),]

write.csv(trainingSet, "Res/trainTestSplit/training_set.csv")
write.csv(testingSet, "Res/trainTestSplit/testing_set.csv")
