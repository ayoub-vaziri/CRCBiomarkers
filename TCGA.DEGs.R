# load packages
library(magrittr)
library(DESeq2)
library(apeglm)
library(org.Hs.eg.db)
library(annotate)

# set working directory
setwd("C:/Users/ayoub/OneDrive/Desktop/TCGA/")

# load data
countData <- read.table("Data/COADREAD.uncv2.mRNAseq_raw_counts.txt", sep = "\t", header = T, row.names = 1) %>%  round(digits = 0)

#### Separating tumor data from normal ####
TP <- which(substr(colnames(countData), 14, 15) == "01")
NT <- which(substr(colnames(countData), 14, 15) == "11")

samples <- rep("X", ncol(countData))
samples[TP] <- "T"
samples[NT] <- "N"

countData = countData[, c(samples != "X")]

samples = samples[c(samples != "X")]

colData <- data.frame(condition = factor(samples))
colData$condition <- relevel(colData$condition, "N")

#### Finding DEGs using DESeq2 package ####
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
memory.limit(20000)
dds <- dds[keep,] %>% DESeq()
resLFC <- lfcShrink(dds, coef = "condition_T_vs_N", type = "apeglm")

#### Converting ID to symbol  ####
symbol <- function(x) strsplit(x, split = '\\|')[[1]][[1]]
id <- function(x) strsplit(x, split = '\\|')[[1]][[2]]
get.symbol <- function(x) lookUp(id(x), 'org.Hs.eg', 'SYMBOL')[[1]]

symbol.id <- rownames(resLFC)
idx <- list(grep("\\?", symbol.id))[[1]]
gene.symbol <- sapply(1:length(symbol.id), function(x) ifelse(x %in% idx, get.symbol(symbol.id[x]), symbol(symbol.id[x])))

#### TCGA_COADREAD top table ####
degs <- data.frame(Gene.symbol = gene.symbol,
                   adj.P.Val = resLFC$padj,
                   P.Value = resLFC$pvalue,
                   logFC = resLFC$log2FoldChange); row.names(degs) <- NULL

degs <- degs[order(degs$P.Value),]
rownames(degs) <- NULL

write.table(degs, file = "Res/TCGA_COADREAD.top.table.tsv", sep = "\t", row.names = F, quote = F)
