library(data.table)
library(RColorBrewer)
library(ggplot2)
library(sva)
library(preprocessCore)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load data ####
###################
exprColorectal <- as.data.frame(fread("Res/mergeDatasets/mergeDatasets.csv")) 
rownames(exprColorectal) <- exprColorectal$V1
exprColorectal <- exprColorectal[,-1]

batchColorectal <- exprColorectal$batch

groupColorectal <- exprColorectal$group

batchGroup <- exprColorectal[,c(1,2)]
batchGroup <- cbind(sample=rownames(exprColorectal), batchGroup)

exprColorectal <- data.matrix(exprColorectal[,-c(1,2)])
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

clrs <- c(
          rep("#1B9E77",length(which(batchColorectal == 1))),
          rep("#666666",length(which(batchColorectal == 2))), 
          rep("#66A61E",length(which(batchColorectal == 3))), 
          rep("#7570B3",length(which(batchColorectal == 4))),
          rep("#A6761D",length(which(batchColorectal == 5)))
          # rep("#D95F02",length(which(batchColorectal == 6)))
          )


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


#### PCA before correction ####
###############################
pcaBeforeCorrection <- prcomp(exprColorectal, scale=FALSE, center = TRUE)

pca.var.per.before <- pcaVarPer(pcaBeforeCorrection)

pc1 <- pcaBeforeCorrection$x[,1]
pc2 <- pcaBeforeCorrection$x[,2]

sdat <- data.frame(x=pc1, y=pc2, batch = grps[batchColorectal])
tdat <- data.frame(sdat[,1:2], grps[batchColorectal], type[groupColorectal])

colnames(tdat)<-c("PC-1","PC-2","Batch","Group")
rownames(tdat)<-rownames(exprColorectal)

write.table(tdat,paste0("Res/batchCorrection/pca_Before_Correction.txt"),quote=TRUE,sep="\t",row.names=TRUE)

png("Res/batchCorrection/pca_Before_Correction.png", width = 2600, height = 2000, res = 300)
fig(sdat, groupColorectal, clrs, pca.var.per.before)
dev.off()
###############################

#### Correct batch by ComBat ####
#################################
combatres <- ComBat(
  dat=t(exprColorectal),
  batch=batchColorectal,
  mod=NULL, 
  par.prior=TRUE, 
  prior.plots=FALSE
)

combatres_ <- as.data.frame(cbind(sample=rownames(t(combatres)), t(combatres)))
combatres_ <- merge(batchGroup, combatres_, by="sample")

write.csv(combatres_, "Res/batchCorrection/corrected_exprColorectal.csv")
write.table(combatres_, paste0("Res/batchCorrection/corrected_exprColorectal.txt"),quote=TRUE,sep="\t",row.names=TRUE)
#################################

#### PCA on combat batch corrected ####
#######################################
pcaAfterCorrection <- prcomp(t(combatres), scale=FALSE, center = TRUE)

pca.var.per.after <- pcaVarPer(pcaAfterCorrection)

pc1 <- pcaAfterCorrection$x[,1]
pc2 <- pcaAfterCorrection$x[,2]

sdat <- data.frame(x=pc1, y=pc2, batch = grps[batchColorectal])
tdat <- data.frame(sdat[,1:2], grps[batchColorectal], type[groupColorectal])

colnames(tdat)<-c("PC-1","PC-2","Batch","Group")
rownames(tdat)<-rownames(exprColorectal)

write.table(tdat,paste0("Res/batchCorrection/pca_After_Correction.txt"),quote=TRUE,sep="\t",row.names=TRUE)

png("Res/batchCorrection/pca_After_Correction.png", width = 2600, height = 2000, res = 300)
fig(sdat, groupColorectal, clrs, pca.var.per.after)
dev.off()
#######################################