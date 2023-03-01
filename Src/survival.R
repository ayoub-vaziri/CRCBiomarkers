library(data.table)
library(survival)
library(survminer)
library(finalfit)
library(forestplot)
library(dplyr)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load and process data ####
###############################
mrnaData <- as.data.frame(fread("Res/survival/Data/coadread_tcga/data_mrna_seq_v2_rsem.txt"))
clinData <- as.data.frame(fread("Res/survival/Data/coadread_tcga/data_clinical_patient.txt"))

mrnaData <- na.omit(mrnaData)
mrnaData <- mrnaData[which(mrnaData$Hugo_Symbol != ""),]
mrnaData <- mrnaData[!duplicated(mrnaData$Hugo_Symbol),]

rownames(mrnaData) <- mrnaData$Hugo_Symbol
mrnaData <- mrnaData[,-c(1,2)]

mrnaData <- log2(mrnaData+1)

colnames(mrnaData) <- gsub("-01", "", colnames(mrnaData))
mrnaData <- as.data.frame(t(mrnaData))
mrnaData <- cbind(PATIENT_ID=rownames(mrnaData), mrnaData)

clinData <- clinData[-c(1:3),]
colnames(clinData) <- clinData[1,]
clinData <- clinData[-1,]
#clinData <- clinData[,c("PATIENT_ID", "OS_STATUS", "OS_MONTHS")]
clinData <- clinData[,c("PATIENT_ID", "DFS_STATUS", "DFS_MONTHS")]
#clinData$OS_STATUS <- gsub(":.*", "", clinData$OS_STATUS)
clinData$DFS_STATUS <- gsub(":.*", "", clinData$DFS_STATUS)
for(x in c(2,3)) clinData[,x] <- as.numeric(clinData[,x])
clinData <- na.omit(clinData)
rownames(clinData) <- NULL

colnames(mrnaData)[4476] <- "CXCL1"
mrnaClin <- merge(clinData, mrnaData, by="PATIENT_ID")
rownames(mrnaClin) <- mrnaClin$PATIENT_ID
mrnaClin <- mrnaClin[,-1]

#write.csv(mrnaClin, "Res/supplementary/OS_mRNA_colorectal.csv")
write.csv(mrnaClin, "Res/supplementary/DFS_mRNA_colorectal.csv")
###############################

#### Survival plot ####
#######################
survplt <- function(fit, dat, gene, title) {
  ggsurvplot(fit = fit, 
             data = dat,
             title = title,
             font.title = c(14, "black"),
             ggtheme = theme_test() + theme(plot.title = element_text(hjust = 0.5)),
             ####### Censor Details ########
             censor = FALSE,
             censor.shape="+",
             censor.size = 6,
             ####### Confidence Intervals ########
             conf.int = FALSE,
             ####### Format Axes #######
             xlab="Months",
             ylab = "Survival probability",
             font.x=c(12), 
             font.y=c(12), 
             font.xtickslab=c(12,"plain"),
             font.ytickslab=c(12,"plain"),
             ######## Format Legend #######
             legend.title = gene,
             font.legend =c(18, "black"),
             legend.labs = c("high", "low"),
             legend = c(.83, .85),
             ######## Plot Dimensions #######
             surv.plot.height = 2,
             surv.plot.width = 2,
             ######## Risk Table #######
             risk.table = FALSE,
             risk.table.height = 0.25,
             ######## p-value details #######
             pval = TRUE,
             pval.method = TRUE,
             pval.size = 7,
             pval.method.size = 5,
             pval.coord = c(1,0.1),
             pval.method.coord = c(1, 0.18)
  )
}
#######################

#### Survival analysis ####
###########################
os_mrna_colorectal <- as.data.frame(fread("Res/supplementary/OS_mRNA_colorectal.csv"))
rownames(os_mrna_colorectal) <- os_mrna_colorectal$V1
os_mrna_colorectal <- os_mrna_colorectal[,-1]

dfs_mrna_colorectal <- as.data.frame(fread("Res/supplementary/DFS_mRNA_colorectal.csv"))
rownames(dfs_mrna_colorectal) <- dfs_mrna_colorectal$V1
dfs_mrna_colorectal <- dfs_mrna_colorectal[,-1]

survAnalysis <- function(gene, dat, time, status, title) {
  data <- dat[, c(status, time, gene)]
  high <- which(data[,3] >= median(data[,3]))
  low <- which(data[,3] < median(data[,3]))
  data[high,3] <- "high"
  data[low,3] <- "low"
  sfit <- survfit(Surv(data[,2], data[,1]) ~ data[,3], data=data)
  return(survplt(sfit, data, gene, title))
}

lassoRfGenes <- fread("Res/venn/lassoRfGenes.txt", header = FALSE)$V1
lassoRfGenes <- lassoRfGenes[-1]
lassoRfGenes <- sort(lassoRfGenes)

Map(function(g) {
  survAnalysis(g, os_mrna_colorectal, "OS_MONTHS", "OS_STATUS", "Overall Survival")
}, lassoRfGenes) -> list_os_plots

Map(function(g) {
  survAnalysis(g, dfs_mrna_colorectal, "DFS_MONTHS", "DFS_STATUS", "Disease Free Survival")
}, lassoRfGenes) -> list_dfs_plots

png(filename = "Res/supplementary/OS_survival_analysis.png", width = 4200, height = 2200, res = 300)
arrange_ggsurvplots(list_os_plots, print = TRUE, nrow = 2, ncol = 4)
dev.off()

png(filename = "Res/supplementary/DFS_survival_analysis.png", width = 4200, height = 2200, res = 300)
arrange_ggsurvplots(list_dfs_plots, print = TRUE, nrow = 2, ncol = 4)
dev.off()

Map(function(g) {
  survAnalysis(g, os_mrna_colorectal, "OS_MONTHS", "OS_STATUS", "Overall Survival")
}, c("BGN", "SPP1")) -> list_os_plots

Map(function(g) {
  survAnalysis(g, dfs_mrna_colorectal, "DFS_MONTHS", "DFS_STATUS", "Disease Free Survival")
}, c("BGN", "SPP1")) -> list_dfs_plots

png(filename = "Res/survival/Res/OS_survival_analysis.png", width = 2600, height = 1200, res = 300)
arrange_ggsurvplots(list_os_plots, print = TRUE, nrow = 1, ncol = 2)
dev.off()

png(filename = "Res/survival/Res/DFS_survival_analysis.png", width = 2600, height = 1200, res = 300)
arrange_ggsurvplots(list_dfs_plots, print = TRUE, nrow = 1, ncol = 2)
dev.off()

list_plots <- c(list_os_plots[1], list_dfs_plots[1], list_os_plots[2], list_dfs_plots[2])

png(filename = "Res/survival/Res/OS_DFS_survival_analysis.png", width = 2800, height = 2600, res = 300)
arrange_ggsurvplots(list_plots, print = TRUE, nrow = 2, ncol = 2)
dev.off()
###########################

#### Univariate Cox regression analysis ####
############################################
explanatory <- lassoRfGenes
dependent_os <- "Surv(OS_MONTHS, OS_STATUS)"

output <- os_mrna_colorectal %>% finalfit(dependent_os, explanatory, add_dependent_label = TRUE)

uniOutput <- data.frame(
  mean = as.numeric(substr(output$`HR (univariable)`, 1, 4)),
  lower = as.numeric(substr(output$`HR (univariable)`, 7, 10)),
  upper = as.numeric(substr(output$`HR (univariable)`, 12, 15)),
  gene = output$`Dependent: Surv(OS_MONTHS, OS_STATUS)`,
  pvalue = substr(output$`HR (univariable)`, 20, 24),
  HR = gsub(",", ")", substr(output$`HR (univariable)`, 1, 16))
)

forestPlot <- uniOutput |> 
  forestplot(labeltext = c(gene, pvalue, HR), 
             #clip = c(-Inf, Inf),
             xticks = c(-0.32, 0, 0.425),
             lwd.xaxis = 1.5,
             lwd.ci = 1.4,
             lwd.zero = 2.2,
             boxsize = 0.2,
             colgap = unit(1, "cm"),
             hrzl_lines = TRUE,
             xlab = "Hazard ratio",
             #title = "Univariate Cox regression",
             vertices = TRUE,
             xlog = TRUE) |> 
  fp_set_style(box = "royalblue4",
               line = "royalblue4",
               align = "lrrr",
               
               hrz_lines = "#999999") |> 
  fp_add_header(gene = c(""),
                pvalue = c("p-value"),
                HR = c("HR (95% CI)"))

png("Res/survival/Res/univariateAnalysis.png", height = 1800, width = 2000, res = 300)
forestPlot
dev.off()
############################################