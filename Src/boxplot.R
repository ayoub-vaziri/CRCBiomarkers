library(dplyr)
library(ggpubr)
library(data.table)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Boxplot for validation data ####
#####################################
train_data <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
bio_data <- fread("Res/exprBatchBioData/bioColorectal.txt", header = FALSE)$V1
lassoGenes <- fread("Res/LASSO/lassoGenes.txt", header = FALSE)$V1

smps <- 1:52

valid_data <- train_data[smps,]
rownames(valid_data) <- valid_data$V1
valid_data <- valid_data[,-1]
valid_data <- valid_data[,lassoGenes]

bio_data <- bio_data[smps]
bio_data <- as.factor(bio_data)
levels(bio_data) <- c("Normal", "Tumor")

valid_data <- cbind(status=bio_data, valid_data)

# pdf boxplots
for(gene in lassoGenes) {
  pdf(paste0("Res/boxplot/pdf/", gene, "_boxplot.pdf"))
  print(
    ggboxplot(valid_data, 
              x = "status",
              y = gene,
              combine = TRUE,
              color = "status", 
              palette = c("red", "blue"),
              ylab = paste(gene, "expression"), 
              xlab = "",
              add = "jitter",                            
              add.params = list(size = 2.5, jitter = 0.2)) + 
      stat_compare_means(comparisons =  list(c("Tumor", "Normal")))
  )
  dev.off()      
}

# png boxplots
cols <- setdiff(colnames(valid_data), "status")
y_labels <- cols

Map(function(x, y) {
  ggboxplot(data = valid_data, x = "status", y = x,
            color = "status",
            palette = c("blue", "red"),
            ylab = paste(y, "expression"),
            xlab = "",
            add = "jitter",
            add.params = list(size = 2, jitter = 0.2)) +
    stat_compare_means(comparisons = list(c("Tumor", "Normal")), size = 5) +
    theme(text = element_text(size = 15),
          legend.text = element_text(size=15, face="bold"))
}, cols, y_labels) -> list_plots

png(filename = "Res/boxplot/boxplot_test.png", width = 4000, height = 4000, res = 300)
ggarrange(plotlist = list_plots, common.legend = TRUE) 
dev.off()
#####################################