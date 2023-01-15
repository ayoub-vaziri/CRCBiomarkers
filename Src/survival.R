library(survival)
library(survminer)
library(data.table)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load and process data ####
###############################
mrnaData <- as.data.frame(fread("Res/survival/Data/coadread_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt"))

mrnaData <- na.omit(mrnaData)
mrnaData <- mrnaData[which(mrnaData$Hugo_Symbol != ""),]
mrnaData <- mrnaData[!duplicated(mrnaData$Hugo_Symbol),]

rownames(mrnaData) <- mrnaData$Hugo_Symbol
mrnaData <- mrnaData[,-c(1,2)]

colnames(mrnaData) <- gsub("-01", "", colnames(mrnaData))
mrnaData <- as.data.frame(t(mrnaData))
mrnaData <- cbind(PATIENT_ID=rownames(mrnaData), mrnaData)

clinData <- as.data.frame(fread("Res/survival/Data/coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt"))

clinData <- clinData[-c(1:3),]
colnames(clinData) <- clinData[1,]
clinData <- clinData[-1,]
clinData <- clinData[,c("PATIENT_ID", "OS_STATUS", "OS_MONTHS")]
clinData <- clinData[-which(clinData$OS_STATUS == ""),]
clinData$OS_STATUS <- gsub(":.*", "", clinData$OS_STATUS)
for(x in c(2,3)) clinData[,x] <- as.numeric(clinData[,x])
rownames(clinData) <- NULL

mrnaClin <- merge(clinData, mrnaData, by="PATIENT_ID")
rownames(mrnaClin) <- mrnaClin$PATIENT_ID
mrnaClin <- mrnaClin[,-1]

write.csv(mrnaClin, "Res/survival/Res/OS_mRNA_colorectal_cell.csv")
###############################

#### Survival plot ####
#######################
survplt <- function(fit, dat, gene) {
  ggsurvplot(fit = fit, 
             data = dat,
             title = "Overall Survival",
             font.title = c(15, "black"),
             ggtheme = theme_test() + theme(plot.title = element_text(hjust = 0.5)),
             ####### Censor Details ########
             censor = FALSE,
             censor.shape="+",
             censor.size = 6,
             ####### Confidence Intervals ########
             conf.int = FALSE,
             ####### Format Axes #######
             xlab="Months",
             ylab = "Percent survival",
             font.x=c(14), 
             font.y=c(14), 
             font.xtickslab=c(14,"plain"),
             font.ytickslab=c(14,"plain"),
             ######## Format Legend #######
             legend.title = gene,
             font.legend =c(13, "black"),
             legend.labs = c("high expression", "low expression"),
             legend = c(.7, .85),
             ######## Plot Dimensions #######
             surv.plot.height = 2,
             surv.plot.width = 2,
             ######## Risk Table #######
             risk.table = FALSE,
             risk.table.height = 0.25,
             ######## p-value details #######
             pval = TRUE,
             pval.method = TRUE,
             pval.size = 5.5,
             pval.method.size = 5,
             pval.coord = c(1,0.1),
             pval.method.coord = c(1, 0.18)
  )
}
#######################

#### Survival analysis ####
###########################
os_mrna_colorectal <- as.data.frame(fread("Res/survival/Res/OS_mRNA_colorectal_cell.csv"))
rownames(os_mrna_colorectal) <- os_mrna_colorectal$V1
os_mrna_colorectal <- os_mrna_colorectal[,-1]

survAnalysis <- function(gene) {
  data <- os_mrna_colorectal[, c("OS_STATUS", "OS_MONTHS", gene)]
  high <- which(data[,3] >= mean(data[,3]))
  low <- which(data[,3] < mean(data[,3]))
  data[high,3] <- "high"
  data[low,3] <- "low"
  sfit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ data[,3], data = data)
  return(survplt(fit = sfit, dat = data, gene = gene))
}

lassoGenes <- c("MYC", "CCND1", "BMP2", "MMP3", "CXCL1", "BGN", "SPP1", "MMP7", "PLAU")

Map(function(g) {
survAnalysis(g)
}, lassoGenes) -> list_plots

png(filename = "Res/survival/Res/OS_survival_analysis.png", width = 4000, height = 4000, res = 300)
arrange_ggsurvplots(list_plots, print = TRUE, nrow = 3, ncol = 3)
dev.off()
###########################