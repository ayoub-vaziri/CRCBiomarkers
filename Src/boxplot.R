library(dplyr)
library(ggpubr)
library(data.table)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Boxplot for validation data ####
#####################################
train_data <- as.data.frame(fread("Res/batchCorrection/corrected_exprColorectal.csv"))
bio_data <- fread("Res/exprBatchBioData/bioColorectal.txt", header = FALSE)$V1

lassoRfGenes <- fread("Res/venn/lassoRfGenes.txt", header = FALSE)$V1
lassoRfGenes <- sort(lassoRfGenes)

smps <- 1:102

valid_data <- train_data[smps,]
rownames(valid_data) <- valid_data$V1
valid_data <- valid_data[,-1]
valid_data <- valid_data[,lassoRfGenes]

bio_data <- bio_data[smps]
bio_data <- as.factor(bio_data)
levels(bio_data) <- c("Normal", "Tumor")

valid_data <- cbind(status=bio_data, valid_data)

# png boxplots
cols <- setdiff(colnames(valid_data), "status")
y_labels <- cols

Map(function(x, y, z) {
  ggboxplot(data = valid_data, x = "status", y = x,
            color = "status",
            palette = c("blue", "red"),
            ylab = paste(y, "expression"),
            xlab = "",
            title = z,
            add = "jitter",
            add.params = list(size = 2, jitter = 0.2)) +
    stat_compare_means(comparisons = list(c("Tumor", "Normal")), size = 5) +
    theme(text = element_text(size = 15),
          title = element_text(size = 18),
          legend.text = element_text(size=15))
}, cols, y_labels, LETTERS[1:9]) -> list_plots

png(filename = "Res/boxplot/boxplot_test.png", width = 5200, height = 3000, res = 300)
ggarrange(plotlist = list_plots, common.legend = TRUE, ncol = 5, nrow = 2) 
dev.off()
#####################################