library(data.table)

# Set the current working directory to the project path
setwd(project_path)

samples <- as.data.frame(fread("Results/DataProcessing/batchCorrection/corrected_exprColorectal.csv"))
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

write.csv(trainingSet, "Results/DataProcessing/trainTestSplit/training_set.csv")
write.csv(testingSet, "Results/DataProcessing/trainTestSplit/testing_set.csv")
