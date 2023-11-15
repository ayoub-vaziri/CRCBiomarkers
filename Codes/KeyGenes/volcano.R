library(limma)
library(data.table)
library(ggrepel)
library(dplyr)
library(readr)

# Set the current working directory to the project path
setwd(project_path)

#### Volcano plot for validation set ####
#########################################
degs <- function(expr, group) {
  design.mat <- model.matrix(~ group + 0)
  colnames(design.mat) <- levels(group)
  lmfit <- lmFit(expr, design.mat)
  cont.mat <- makeContrasts(contrasts = "Tumor-Normal", levels = design.mat)
  cont.fit <- contrasts.fit(lmfit, cont.mat)
  fit <- eBayes(cont.fit, 0.01)
  toptable <- topTable(fit = fit, number = Inf, adjust.method = "fdr", sort.by="B")
  return(toptable)
}

expr <- as.data.frame(fread("Results/DataProcessing/trainTestSplit/training_set.csv"))
bio <- expr$group

rownames(expr) <- expr$V1
expr <- expr[,-c(1,2)]
expr <- as.data.frame(t(expr))

bio <- gsub(1, "N", bio)
bio <- gsub(2, "C", bio)
bio <- as.factor(bio)
levels(bio) <- c("Tumor", "Normal")

exprData <- expr
bioData <- bio

tT <- degs(exprData, bioData)
tT$Gene.symbol <- rownames(tT)

data <- tT %>% 
  mutate(
    Regulated = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "Up",
                    logFC < -1 & adj.P.Val < 0.05 ~ "Down",
                    TRUE ~ "Not")
  ) 

p <- ggplot(data, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color = Regulated), size = 1) +
  xlab(expression("log2(FoldChange)")) + 
  ylab(expression("-log10(adj.P.Val)")) +
  scale_color_manual(values = c("blue", "black", "red")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

data %>% 
  count(Regulated)

top_genes <- data[data$adj.P.Val < 1e-20 & abs(data$logFC) > 1.4, ]

options(ggrepel.max.overlaps = 15)

p <- p +
  geom_label_repel(data = top_genes,
                   mapping = aes(logFC, -log10(adj.P.Val), label = Gene.symbol),
                   size = 2) +
  theme_bw()

theme <- theme(strip.text.y = element_text(),
               axis.text = element_text(colour = "black", size=12),
               axis.title.x = element_text(colour = "black", face="bold", size=14),
               axis.title.y = element_text(colour = "black", face="bold", size=14),
               axis.text.x = element_text(colour = "black", size=12),
               legend.background = element_rect(fill = "white"),
               panel.grid.major = element_line(colour = "white"),
               legend.title=element_text(colour="black",size=12),
               legend.text=element_text(colour="black",size=12),
               strip.text.x = element_text(colour = "black", size = 12),
               plot.title=element_text(size=12, face="bold",hjust=.5),
               panel.border=element_rect(colour="black",size=1.2),
               legend.key=element_blank())

png("Results/KeyGenes/volcano/volcano.png", height = 2000, width = 2400, res = 300)
p + theme
dev.off()
#########################################