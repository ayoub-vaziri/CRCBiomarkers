library(data.table)
library(factoextra)

# Set the current working directory to the project path
setwd(project_path)

#### PCA plot for validation set ####
#####################################
train_data <- as.data.frame(fread("Results/DataProcessing/trainTestSplit/training_set.csv"))
test_data <- as.data.frame(fread("Results/DataProcessing/trainTestSplit/testing_set.csv"))

colnames(train_data)[1] <- "sample"
colnames(test_data)[1] <- "sample"

train_group <- train_data$group
test_group <- test_data$group

train_group <- as.factor(train_group)
levels(train_group) <- c("Normal", "Tumor")

test_group <- as.factor(test_group)
levels(test_group) <- c("Normal", "Tumor")

rownames(train_data) <- train_data$sample
rownames(test_data) <- test_data$sample

train_data <- train_data[,-c(1,2)]
test_data <- test_data[,-c(1,2)]

type <- c("Normal", "Tumor")

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

fig <- function(sdat, group, clrs) {
  fch=c(15,16)
  Group=factor(fch[group],labels=type)
  ggplot(data=sdat, aes(x, y,colour=Group)) +
    geom_point(aes(shape=Group),size=5.5,alpha=.7) +
    geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.6) +
    geom_vline(xintercept=0, linetype="dashed", colour= "black", size = 0.6) +
    labs(title="",x="PC-1",y="PC-2") +
    scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=5))) +
    scale_colour_manual(labels=type,values=c("cyan3", "coral2")) +
    theme_bw() +
    xlab("PC1") +
    ylab("PC2") +
    theme
}

train_clrs <- c(
  rep("cyan3",length(which(train_group == "Normal"))),
  rep("coral2",length(which(train_group == "Tumor")))
)

test_clrs <- c(
  rep("cyan3",length(which(test_group == "Normal"))),
  rep("coral2",length(which(test_group == "Tumor")))
)

train_pca <- prcomp(train_data, scale = F, center = T)

pc1 <- train_pca$x[,1]
pc2 <- train_pca$x[,2]

sdat <- data.frame(x=pc1, y=pc2, group=type[as.numeric(train_group)])

png("Results/DataProcessing/PCA/train.pcaplot.png", width = 2600, height = 2000, res = 300)
fig(sdat, as.numeric(train_group), train_clrs)
dev.off()


test_pca <- prcomp(test_data, scale = F, center = T)

pc1 <- test_pca$x[,1]
pc2 <- test_pca$x[,2]

sdat <- data.frame(x=pc1, y=pc2, group=type[as.numeric(test_group)])

png("Results/DataProcessing/PCA/test.pcaplot.png", width = 2600, height = 2000, res = 300)
fig(sdat, as.numeric(test_group), train_clrs)
dev.off()