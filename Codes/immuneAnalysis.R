library(IOBR)
library(ggplot2)
library(corrplot)

setwd(project_path)

train <- as.data.frame(data.table::fread("Results/batchCorrection/corrected_training_set.csv"))
valid <- as.data.frame(data.table::fread("Results/batchCorrection/corrected_validation_set.csv"))

colnames(train)[c(2,4)] <- c("ID", "Group")
colnames(valid)[c(2,4)] <- c("ID", "Group")

train <- train[,-c(1,3)]
valid <- valid[,-c(1,3)] 

dat <- rbind(train, valid)
eset <- dat
rownames(eset) <- eset$ID
eset <- eset[,-c(1,2)]

cibersort <- deconvo_tme(eset = as.data.frame(t(eset)), 
                         method = "cibersort", 
                         arrays = TRUE, 
                         perm = 1000)

#### violin plot ####
#####################
cibersort <- merge(dat[,c(1,2)], cibersort, by = "ID")
cibersort$Group <- gsub(1, "Non", cibersort$Group)
cibersort$Group <- gsub(2, "CRC", cibersort$Group)
colnames(cibersort) <- gsub("_CIBERSORT", "", colnames(cibersort))
cibersort$Group <- as.factor(cibersort$Group)

cibersort_long <- reshape2::melt(cibersort[,-c(1, 25:27)], 
                                 id.vars = "Group", 
                                 variable.name = "Cell", 
                                 value.name = "Fraction")

p_values <- sapply(levels(cibersort_long$Cell), function(cell) {
  subset_data <- subset(cibersort_long, Cell == cell)
  wilcox.test(Fraction ~ Group, data = subset_data)$p.value
})

maxFrac <- sapply(levels(cibersort_long$Cell), function(cell) {
  subset_data <- subset(cibersort_long, Cell == cell)
  max(subset_data$Fraction) + 0.02
})

annotations <- data.frame(
  Cell = names(p_values),
  p_value = p_values,
  y = maxFrac
)

p <- ggplot(cibersort_long, aes(x = Cell, y = Fraction, fill = Group)) + 
  geom_violin(scale = "width", adjust = 1) + 
  scale_fill_manual(values = c("red", "blue")) +
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(0.8), color="grey") +
  labs(x = "", y = "Fraction", fill = "Group") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        panel.border=element_rect(colour="black",size=1.5),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
        )

png(filename = "Results/immuneAnalysis/vioplot.png", width = 3800, height = 2400, res = 300)
p + annotate("text",  
             x = annotations$Cell, y = annotations$y, 
             label = ifelse(annotations$p_value < 0.001, 
                            "p<0.001", 
                            paste0("p=", round(annotations$p_value, 3))),
             hjust = 0.5, vjust = 0.5, size = 3.6) +
  annotate("text", 
           x = annotations$Cell, y = annotations$y-0.019, 
           label = "───", size = 2,
           hjust = 0.5, vjust = 0.5)
dev.off()
#####################

#### Lollipop plot ####
#######################
biomarkers <- data.table::fread("Results/LASSO/lassoGenes.txt", header = FALSE)$V1
biomarkers <- sort(biomarkers)
biomarkers <- biomarkers[-c(3,5)]

biomarkers_df <- subset(dat, select = c("ID", biomarkers))
cibersort_df <- cibersort[,c(1, 3:24)]
df <- merge(biomarkers_df, cibersort_df, by = "ID")

correlations <- lapply(biomarkers, function(marker) {
  cor_results <- lapply(colnames(df)[11:32], function(cell) {
    cor_test <- cor.test(df[, marker], df[, cell], method = "spearman")
    cor_df <- data.frame(Cell = cell, Correlation = cor_test$estimate, p_value = cor_test$p.value)
    cor_df$Gene <- marker
    return(cor_df)
  })
  return(do.call(rbind, cor_results))
})

i = 9

correlation <- correlations[[i]]
gene <- unique(correlation$Gene)

corr <- correlation[order(correlation$Correlation, decreasing = T),]
rownames(corr) <- NULL

corr$Cell <- gsub("_", " ", corr$Cell)
corr$Cell <- forcats::fct_reorder(corr$Cell, corr$Correlation, .desc = FALSE)

corr$p_value_category <- cut(corr$p_value, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = FALSE)

format_p_value <- function(p_value) {
  ifelse(p_value < 0.001, "<0.001", round(p_value, 3))
}

png(filename = paste0("Results/immuneAnalysis/", gene, ".png"), width = 2400, height = 2000, res = 300)
ggplot(corr, aes(x = Cell, y = Correlation, size = abs(Correlation), fill = p_value)) +
  geom_segment(aes(x = Cell, xend = Cell, y = 0, yend = Correlation), color = "black", size = 1.2) +
  geom_point(shape = 21, color = "black") +
  scale_size_continuous(range = c(3, 8), limits = c(0, 1), breaks = seq(0.1, 0.5, by = 0.1)) +
  scale_fill_continuous(limits = c(0, 1),
                    labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"),
                    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 1), breaks = seq(-0.8, 0.8, by = 0.4)) +
  geom_text(aes(label = format_p_value(p_value)),
            x = as.numeric(corr$Cell), y = Inf,
            color = ifelse(corr$p_value < 0.05, "red", "black"),
            hjust = 0, size = 5) +
  coord_flip(clip = "off") +
  labs(fill = "p-value", size = "Correlation") +
  labs(x = "", y = "Correlation Coefficient") +
  labs(size = "aps(cor)") +
  labs(fill = "pvalue") +
  labs(title = gene) +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 17, color = "black"),
    legend.margin = margin(0, 40, 0, 0),
    legend.box.margin = margin(0, 0, 0, 40),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),

  ) +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "black")))
dev.off()
#######################

#### Heatmap ####
#################
df <- cibersort[, c(3:24)]
colnames(df) <- gsub("_", " ", colnames(df))
cor_matrix <- cor(df, method = "pearson")

color_palette <- colorRampPalette(c("red", "white", "blue"))(200)

png(filename = paste0("Results/immuneAnalysis/heatmap.png"), width = 2200, height = 2200, res = 300)
corrplot(cor_matrix, 
         method = "color", 
         addrect = 2, 
         addCoef.col = TRUE, 
         tl.cex = 0.8,
         tl.col = "black",
         number.cex = 0.5,
         type = "full",
         col = rev(color_palette)
)
dev.off()
#################