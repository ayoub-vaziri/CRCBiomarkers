library(data.table)
library(RColorBrewer)
library(ggplot2)
library(sva)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load data ####
###################
exprColorectal <- as.data.frame(fread("Res/exprBatchBioData/exprColorectal.csv")) 
rownames(exprColorectal) <- exprColorectal$V1
exprColorectal <- exprColorectal[,-1]

exprColorectal <- data.matrix(exprColorectal)

batchColorectal <- as.numeric(fread("Res/exprBatchBioData/batchColorectal.txt")$V1)

bioColorectal <- as.numeric(fread("Res/exprBatchBioData/bioColorectal.txt")$V1)
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

clrs <- c(rep("#1B9E77",nrow(expr25070)), 
          rep("#66A61E",nrow(expr25071)), 
          rep("#7570B3",nrow(expr44076)), 
          rep("#D95F02",nrow(expr106582)),
          rep("#E6AB02",nrow(expr110223)),
          rep("#E7298A",nrow(expr110224))
          )

grps <- c("GSE25070", "GSE25071", "GSE44076", "GSE106582", "GSE110223", "GSE110224")

type <- c("Normal", "Tumor")

values <- sort(brewer.pal(n=8,name="Dark2")[1:6])
########################

#### PC figure ####
###################
fig <- function(sdat) {
  fch=c(15,16)
  Group=factor(fch[bioColorectal],labels=type)
  Batch=clrs
  ggplot(data=sdat, aes(x, y,colour=Batch)) +
    geom_point(aes(shape=Group),size=5.5,alpha=.7) +
    labs(title="",x="PC-1",y="PC-2") +
    scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=5))) +
    scale_colour_manual(labels=grps,values=values) +
    theme_bw() +
    theme
}
###################

#### PCA before correction ####
###############################
pcaBeforeCorrection <- prcomp(exprColorectal,scale=TRUE)

x <- pcaBeforeCorrection$x[,1]
y <- pcaBeforeCorrection$x[,2]

sdat <- data.frame(x=x, y=y, group = grps[batchColorectal])
tdat <- data.frame(sdat[,1:2], grps[batchColorectal], type[bioColorectal])

colnames(tdat)<-c("PC-1","PC-2","Batch","Group")
rownames(tdat)<-rownames(exprColorectal)

write.table(tdat,paste0("Res/batchCorrection/pca_Before_Correction.txt"),quote=TRUE,sep="\t",row.names=TRUE)

png("Res/batchCorrection/pca_Before_Correction.png", width = 2600, height = 2000, res = 300)
fig(sdat)
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

write.csv(t(combatres), "Res/batchCorrection/corrected_exprColorectal.csv")
write.table(t(combatres),paste0("Res/batchCorrection/corrected_exprColorectal.txt"),quote=TRUE,sep="\t",row.names=TRUE)
#################################

#### PCA on combat batch corrected ####
#######################################
pcaAfterCorrection <- prcomp(t(combatres),scale=TRUE)

x <- pcaAfterCorrection$x[,1]
y <- pcaAfterCorrection$x[,2]

sdat <- data.frame(x=x, y=y, group = grps[batchColorectal])
tdat <- data.frame(sdat[,1:2], grps[batchColorectal], type[bioColorectal])

colnames(tdat)<-c("PC-1","PC-2","Batch","Group")
rownames(tdat)<-rownames(exprColorectal)

write.table(tdat,paste0("Res/batchCorrection/pca_After_Correction.txt"),quote=TRUE,sep="\t",row.names=TRUE)

png("Res/batchCorrection/pca_After_Correction.png", width = 2600, height = 2000, res = 300)
fig(sdat)
dev.off()
#######################################