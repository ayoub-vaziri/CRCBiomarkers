library(data.table)
library(RColorBrewer)
library(ggplot2)
library(sva)
library(preprocessCore)

# Set the current working directory to the project path
setwd(project_path)

#### Load data ####
###################
training_set <- as.data.frame(fread("Results/trainValidSplit/training_set.csv"))
validation_set <- as.data.frame(fread("Results/trainValidSplit/validation_set.csv"))

rownames(training_set) <- training_set$V1
rownames(validation_set) <- validation_set$V1

training_set <- training_set[,-1]
validation_set <- validation_set[,-1]

train_batch <- training_set$batch
valid_batch <- validation_set$batch

train_group <- training_set$group
valid_group <- validation_set$group

train_batch_group <- training_set[,c(1,2)]
valid_batch_group <- validation_set[,c(1,2)]

train_batch_group <- cbind(sample=rownames(training_set), train_batch_group)
valid_batch_group <- cbind(sample=rownames(validation_set), valid_batch_group)

training_data <- data.matrix(training_set[,-c(1,2)])
validation_data <- data.matrix(validation_set[,-c(1,2)])

###################

#### Set parameters ####
########################
theme <- theme(strip.text.y = element_text(),
               axis.text = element_text(colour = "black", size=15),
               axis.title.x = element_text(colour = "black", face="bold", size=15),
               axis.title.y = element_text(colour = "black", face="bold", size=15),
               axis.text.x = element_text(colour = "black", size=15),
               legend.background = element_rect(fill = "white"),
               panel.grid.major = element_line(colour = "white"),
               legend.title=element_text(colour="black",size=15),
               legend.text=element_text(colour="black",size=15),
               strip.text.x = element_text(colour = "black", size = 15),
               plot.title=element_text(size=15, face="bold",hjust=.5,margin=margin(b=2,unit="pt")),
               panel.border=element_rect(colour="black",size=2),
               legend.key=element_blank())

values <- sort(brewer.pal(n=5,name="Dark2")[1:5])

grps <- c("GSE10950", "GSE25070", "GSE41328", "GSE74602", "GSE142279")

type <- c("Normal", "Tumor")
########################

#### PC figure ####
###################
fig <- function(sdat, group, clrs, pca.var.per) {
  fch=c(15,16)
  Group=factor(fch[group],labels=type)
  Batch=clrs
  ggplot(data=sdat, aes(x, y,colour=Batch)) +
    geom_point(aes(shape=Group),size=5.5,alpha=.7) +
    labs(title="",x="PC-1",y="PC-2") +
    scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=5))) +
    scale_colour_manual(labels=grps,values=values) +
    theme_bw() + 
    xlab(paste0("PC1: ", pca.var.per[1], "% variance")) +
    ylab(paste0("PC2: ", pca.var.per[2], "% variance")) +
    theme
}
###################

#### Calculate the percentage variance ####
###########################################
pcaVarPer <- function(pca) {
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var / sum(pca.var) * 100, 2)
  return(pca.var.per)
}
###########################################


#### PCA before correction for training set ####
################################################
pcaBeforeCorrection <- prcomp(training_data, scale=FALSE, center = TRUE)

pca.var.per.before <- pcaVarPer(pcaBeforeCorrection)

pc1 <- pcaBeforeCorrection$x[,1]
pc2 <- pcaBeforeCorrection$x[,2]

sdat <- data.frame(x=pc1, y=pc2, batch = grps[train_batch])
tdat <- data.frame(sdat[,1:2], grps[train_batch], type[train_group])

colnames(tdat)<-c("PC-1","PC-2","Batch","Group")
rownames(tdat)<-rownames(training_data)

write.table(tdat,paste0("Results/batchCorrection/pca_Before_Correction_for_training_set.txt"),quote=TRUE,sep="\t",row.names=TRUE)

clrs <- c(
  rep("#1B9E77",length(which(train_batch == 1))),
  rep("#666666",length(which(train_batch == 2))), 
  rep("#66A61E",length(which(train_batch == 3))), 
  rep("#7570B3",length(which(train_batch == 4))),
  rep("#A6761D",length(which(train_batch == 5)))
  # rep("#D95F02",length(which(batchColorectal == 6)))
)

png("Results/batchCorrection/pca_Before_Correction_for_training_set.png", width = 2600, height = 2000, res = 300)
fig(sdat, train_group, clrs, pca.var.per.before)
dev.off()
################################################

#### PCA before correction for validation set ####
################################################
pcaBeforeCorrection <- prcomp(validation_data, scale=FALSE, center = TRUE)

pca.var.per.before <- pcaVarPer(pcaBeforeCorrection)

pc1 <- pcaBeforeCorrection$x[,1]
pc2 <- pcaBeforeCorrection$x[,2]

sdat <- data.frame(x=pc1, y=pc2, batch = grps[valid_batch])
tdat <- data.frame(sdat[,1:2], grps[valid_batch], type[valid_group])

colnames(tdat)<-c("PC-1","PC-2","Batch","Group")
rownames(tdat)<-rownames(validation_data)

write.table(tdat,paste0("Results/batchCorrection/pca_Before_Correction_for_validation_set.txt"),quote=TRUE,sep="\t",row.names=TRUE)

clrs <- c(
  rep("#1B9E77",length(which(valid_batch == 1))),
  rep("#666666",length(which(valid_batch == 2))), 
  rep("#66A61E",length(which(valid_batch == 3))), 
  rep("#7570B3",length(which(valid_batch == 4))),
  rep("#A6761D",length(which(valid_batch == 5)))
  # rep("#D95F02",length(which(batchColorectal == 6)))
)

png("Results/batchCorrection/pca_Before_Correction_for_validation_set.png", width = 2600, height = 2000, res = 300)
fig(sdat, valid_group, clrs, pca.var.per.before)
dev.off()
################################################

#### Correct batch by ComBat for training set ####
##################################################
train_combatres <- ComBat(
  dat=t(training_data),
  batch=train_batch,
  mod=NULL, 
  par.prior=TRUE, 
  prior.plots=FALSE
)

combatres_ <- as.data.frame(cbind(sample=rownames(t(train_combatres)), t(train_combatres)))
combatres_ <- merge(train_batch_group, combatres_, by="sample")

write.csv(combatres_, "Results/batchCorrection/corrected_training_set.csv")
write.table(combatres_, paste0("Results/batchCorrection/corrected_training_set.txt"),quote=TRUE,sep="\t",row.names=TRUE)
##################################################

#### Correct batch by ComBat for training set ####
##################################################
valid_combatres <- ComBat(
  dat=t(validation_data),
  batch=valid_batch,
  mod=NULL, 
  par.prior=TRUE, 
  prior.plots=FALSE
)

combatres_ <- as.data.frame(cbind(sample=rownames(t(valid_combatres)), t(valid_combatres)))
combatres_ <- merge(valid_batch_group, combatres_, by="sample")

write.csv(combatres_, "Results/batchCorrection/corrected_validation_set.csv")
write.table(combatres_, paste0("Results/batchCorrection/corrected_validation_set.txt"),quote=TRUE,sep="\t",row.names=TRUE)
##################################################

#### PCA on combat batch corrected for training set ####
########################################################
pcaAfterCorrection <- prcomp(t(train_combatres), scale=FALSE, center = TRUE)

pca.var.per.after <- pcaVarPer(pcaAfterCorrection)

pc1 <- pcaAfterCorrection$x[,1]
pc2 <- pcaAfterCorrection$x[,2]

sdat <- data.frame(x=pc1, y=pc2, batch = grps[train_batch])
tdat <- data.frame(sdat[,1:2], grps[train_batch], type[train_group])

colnames(tdat)<-c("PC-1","PC-2","Batch","Group")
rownames(tdat)<-rownames(training_data)

write.table(tdat,paste0("Results/batchCorrection/pca_After_Correction_for_training_set.txt"),quote=TRUE,sep="\t",row.names=TRUE)

clrs <- c(
  rep("#1B9E77",length(which(train_batch == 1))),
  rep("#666666",length(which(train_batch == 2))), 
  rep("#66A61E",length(which(train_batch == 3))), 
  rep("#7570B3",length(which(train_batch == 4))),
  rep("#A6761D",length(which(train_batch == 5)))
  # rep("#D95F02",length(which(batchColorectal == 6)))
)

png("Results/batchCorrection/pca_After_Correction_training_set.png", width = 2600, height = 2000, res = 300)
fig(sdat, train_group, clrs, pca.var.per.after)
dev.off()
########################################################


#### PCA on combat batch corrected for validation set ####
##########################################################
pcaAfterCorrection <- prcomp(t(valid_combatres), scale=FALSE, center = TRUE)

pca.var.per.after <- pcaVarPer(pcaAfterCorrection)

pc1 <- pcaAfterCorrection$x[,1]
pc2 <- pcaAfterCorrection$x[,2]

sdat <- data.frame(x=pc1, y=pc2, batch = grps[valid_batch])
tdat <- data.frame(sdat[,1:2], grps[valid_batch], type[valid_group])

colnames(tdat)<-c("PC-1","PC-2","Batch","Group")
rownames(tdat)<-rownames(validation_data)

write.table(tdat,paste0("Results/batchCorrection/pca_After_Correction_for_validation_set.txt"),quote=TRUE,sep="\t",row.names=TRUE)

clrs <- c(
  rep("#1B9E77",length(which(valid_batch == 1))),
  rep("#666666",length(which(valid_batch == 2))), 
  rep("#66A61E",length(which(valid_batch == 3))), 
  rep("#7570B3",length(which(valid_batch == 4))),
  rep("#A6761D",length(which(valid_batch == 5)))
  # rep("#D95F02",length(which(batchColorectal == 6)))
)

png("Results/batchCorrection/pca_After_Correction_validation_set.png", width = 2600, height = 2000, res = 300)
fig(sdat, valid_group, clrs, pca.var.per.after)
dev.off()
##########################################################