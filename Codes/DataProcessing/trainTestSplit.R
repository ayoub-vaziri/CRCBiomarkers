library(data.table)

# Set the current working directory to the project path
setwd(project_path)

samples <- as.data.frame(fread("Results/DataProcessing/batchCorrection/corrected_exprColorectal.csv"))
rownames(samples) <- samples$sample
samples <- samples[,-c(1,2,3)]

cancerSamples <- rownames(samples)[which(samples$group == 2)]
normalSamples <- rownames(samples)[which(samples$group == 1)]

trainingCancerSamples <- cancerSamples[sample(1:length(cancerSamples), ceiling(0.7*length(cancerSamples)))]
validationCancerSamples <- cancerSamples[-which(cancerSamples %in% trainingCancerSamples)]

trainingNormalSamples <- normalSamples[sample(1:length(normalSamples), ceiling(0.7*length(normalSamples)))]
validationNormalSamples <- normalSamples[-which(normalSamples %in% trainingNormalSamples)]

trainingSamples <- union(trainingCancerSamples, trainingNormalSamples)
validationSamples <- union(validationCancerSamples, validationNormalSamples)

trainingSet <- samples[which(rownames(samples) %in% trainingSamples),]
validationSet <- samples[which(rownames(samples) %in% validationSamples),]

write.csv(trainingSet, "Results/DataProcessing/trainTestSplit/training_set.csv")
write.csv(validationSet, "Results/DataProcessing/trainTestSplit/validation_set.csv")
