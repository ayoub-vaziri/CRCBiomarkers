library(dplyr)
library(ggpubr)
library(data.table)

# Set the current working directory to the project path
setwd("PROJECT_PATH")

#### Boxplot for test data ####
#####################################
validation <- as.data.frame(fread("Results/DiagnosticGenes/trainTestSplit/validation_set.csv"))

rownames(validation) <- validation$V1
group <- validation$group
group <- as.factor(group)
levels(group) <- c("Normal", "Tumor")
validation <- validation[,-c(1,2)]

lassoGenes <- fread("Results/DiagnosticGenes/LASSO/lassoGenes.txt", header = FALSE)$V1
lassoGenes <- sort(lassoGenes)
lassoGenes <- lassoGenes[-c(3,10)]

validation <- validation[, lassoGenes]
test <- cbind(group=group, validation)

# png boxplots
cols <- setdiff(colnames(test), "group")
y_labels <- cols

Map(function(x, y, z) {
  ggboxplot(data = test, x = "group", y = x,
            color = "group",
            palette = c("blue", "red"),
            ylab = paste(y, "expression"),
            xlab = "",
            title = z,
            add = "jitter",
            add.params = list(size = 2, jitter = 0.2)) +
    stat_compare_means(comparisons = list(c("Tumor", "Normal")), method="wilcox", label = "p.format", size = 5) +
    theme(text = element_text(size = 15),
          title = element_text(size = 18),
          legend.text = element_text(size=15)) 
}, cols, y_labels, LETTERS[1:8]) -> list_plots

png(filename = "Results/DiagnosticGenes/boxplot/boxplot_test.png", width = 4500, height = 3000, res = 300)
ggarrange(plotlist = list_plots, common.legend = TRUE, ncol = 4, nrow = 2) 
dev.off()
#####################################