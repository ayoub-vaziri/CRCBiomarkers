library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(annotate)
library(preprocessCore)
library(EnsDb.Hsapiens.v79)

# Set the current working directory to the project path
setwd(project_path)

#### Load data ####
###################
gse10950 <- getGEO("GSE10950", destdir = "Data/")
gse25070 <- getGEO("GSE25070", destdir = "Data/")
gse41328 <- getGEO("GSE41328", destdir = "Data/")
gse74602 <- getGEO("GSE74602", destdir = "Data/")
gse142279 <- getGEO("GSE142279", destdir = "Data/")

gse10950 <- gse10950[[1]]
gse25070 <- gse25070[[1]]
gse41328 <- gse41328[[1]]
gse74602 <- gse74602[[1]]
gse142279 <- gse142279[[1]]

phen10950 <- pData(gse10950)
phen25070 <- pData(gse25070)
phen41328 <- pData(gse41328)
phen74602 <- pData(gse74602)
phen142279 <- pData(gse142279)

feat10950 <- fData(gse10950)
feat25070 <- fData(gse25070)
feat41328 <- fData(gse41328)
feat74602 <- fData(gse74602)
feat142279 <- fData(gse142279)

# GSE10950 preprocessing
expr10950 <- as.data.frame(fread("Data/GSE10950_Primary_Tumor_nonorm_nobkgd_GEOarchive_Matrix.txt.gz"))
colnames(expr10950) <- expr10950[2,]
expr10950 <- expr10950[-c(1,2),]
rownames(expr10950) <- expr10950$ID_REF
expr10950 <- expr10950[,-1]
for(x in colnames(expr10950)) expr10950[,x] <- as.numeric(expr10950[,x])
expr10950 <- expr10950[,grep("Signal", colnames(expr10950))]
colnames(expr10950) <- gsub(".AVG_Signal", "", colnames(expr10950))
phen10950 <- phen10950[,c(2,8,19)]
colnames(phen10950) <- c("sample", "ID", "group")
phen10950$group <- gsub("Normal", "N", phen10950$group)
phen10950$group <- gsub("Primay tumor", "T", phen10950$group)
expr10950 <- as.data.frame(t(expr10950))
expr10950 <- cbind(ID=rownames(expr10950), expr10950)
expr10950 <- merge(phen10950, expr10950, by="ID")
rownames(expr10950) <- expr10950$sample
expr10950 <- expr10950[,-c(1,2,3)]
expr10950 <- as.data.frame(t(expr10950))
feat10950 <- feat10950[,c(1,12)]
colnames(feat10950) <- c("ID", "symbol")
expr10950 <- log2(expr10950+1)
expr10950 <- cbind(ID=rownames(expr10950), expr10950)
expr10950 <- merge(feat10950, expr10950, by="ID")

# GSE74602 preprocessing
ex74602 <- as.data.frame(fread("Data/GSE74602_non_normalized.txt.gz"))
rownames(ex74602) <- ex74602$ID_REF
ex74602 <- ex74602[,-1]
ex74602 <- as.data.frame(t(ex74602))
ex74602 <- cbind(title=rownames(ex74602), ex74602)
ph74602 <- phen74602[,c(1,2)]
ex74602 <- merge(ph74602, ex74602, by="title")
rownames(ex74602) <- ex74602$geo_accession
ex74602 <- ex74602[,-c(1,2)]
ex74602 <- as.data.frame(t(ex74602))

expr25070 <- exprs(gse25070)
expr41328 <- log2(exprs(gse41328))
expr74602 <- log2(ex74602)

# GSE142279 preprocessing
expr142279 <- as.data.frame(fread("Data/GSE142279_FPKM.xls.gz"))
rownames(expr142279) <- expr142279$ID
expr142279 <- expr142279[,-1]
expr142279 <- log2(expr142279+1)
expr142279 <- as.data.frame(t(expr142279))
expr142279 <- cbind(title=rownames(expr142279), expr142279)
expr142279 <- merge(phen142279[,c(1,2)], expr142279, by="title")
rownames(expr142279) <- expr142279$geo_accession
expr142279 <- as.data.frame(t(expr142279[,-c(1,2)]))
expr142279 <- cbind(ID=rownames(expr142279), expr142279)
symbols142279 <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                   keys= expr142279$ID, 
                                   keytype = "GENEID", 
                                   columns = c("GENEID", "SYMBOL")
                                   )
colnames(symbols142279) <- c("ID", "symbol")
expr142279 <- merge(symbols142279, expr142279, by="ID")
###########

# common genes 
commonGenes <- Reduce(intersect, list(expr10950$symbol,
                                      feat25070$Symbol, 
                                      symbol41328,
                                      feat74602$Symbol,
                                      expr142279$symbol
                                      )
                      )

editSymbol <- function(featData, geneColName) {
  rows <- which(unlist(lapply(strsplit(featData[,geneColName], " /// "), function(x) ifelse(length(x) > 1, TRUE, FALSE))))
  for(x in rows) {
    ss <- strsplit(featData[,geneColName][x], " /// ")[[1]]
    idx <- which(ss %in% commonGenes)[1]
    featData[,geneColName][x] <- ss[idx[1]]
  }
  return(featData)
}

feat25070 <- feat25070[,c(1,12)]
colnames(feat25070) <- c("ID", "symbol")

feat41328 <- editSymbol(feat41328, "Gene Symbol")
feat41328 <- feat41328[,c(1,11)]
colnames(feat41328) <- c("ID", "symbol")

feat74602 <- feat74602[,c(1,12)]
colnames(feat74602) <- c("ID", "symbol")

expr25070 <- cbind(ID=rownames(expr25070), expr25070) %>% as.data.frame()
expr41328 <- cbind(ID=rownames(expr41328), expr41328) %>% as.data.frame()
expr74602 <- cbind(ID=rownames(expr74602), expr74602) %>% as.data.frame()

expr25070 <- merge(feat25070, expr25070, by = "ID")
expr41328 <- merge(feat41328, expr41328, by = "ID")
expr74602 <- merge(feat74602, expr74602, by = "ID")

expr10950 <- expr10950[which(expr10950$symbol != ""),]
expr25070 <- expr25070[which(expr25070$symbol != ""),]
expr41328 <- expr41328[which(expr41328$symbol != ""),]
expr74602 <- expr74602[which(expr74602$symbol != ""),]
expr142279 <- expr142279[which(expr142279$symbol != ""),]

expr10950 <- expr10950[!duplicated(expr10950$symbol),]
expr25070 <- expr25070[!duplicated(expr25070$symbol),]
expr41328 <- expr41328[!duplicated(expr41328$symbol),]
expr74602 <- expr74602[!duplicated(expr74602$symbol),]
expr142279 <- expr142279[!duplicated(expr142279$symbol),]

rownames(expr10950) <- expr10950$symbol
rownames(expr25070) <- expr25070$symbol
rownames(expr41328) <- expr41328$symbol
rownames(expr74602) <- expr74602$symbol
rownames(expr142279) <- expr142279$symbol

expr10950 <- expr10950[,-c(1,2)]
expr25070 <- expr25070[,-c(1,2)]
expr41328 <- expr41328[, -c(1,2)]
expr74602 <- expr74602[,-c(1,2)]
expr142279 <- expr142279[,-c(1,2)]

expr10950 <- na.omit(expr10950)
expr25070 <- na.omit(expr25070)
expr41328 <- na.omit(expr41328)
expr74602 <- na.omit(expr74602)
expr142279 <- na.omit(expr142279)

write.csv(expr10950, "Results/mergeDatasets/expr10950.csv")
write.csv(expr25070, "Results/mergeDatasets/expr25070.csv")
write.csv(expr41328, "Results/mergeDatasets/expr41328.csv")
write.csv(expr74602, "Results/mergeDatasets/expr74602.csv")
write.csv(expr142279, "Results/mergeDatasets/expr142279.csv")
######################

#### Identifying biological classes ####
########################################
for(x in 1:nrow(phen25070)) {
  status <- strsplit(phen25070$title[x], split = " ")[[1]][1]
  if(status == "Colorectal") phen25070$title[x] <- "T" else phen25070$title[x] <- "N"
}

for(x in 1:nrow(phen41328)) {
  status <- strsplit(phen41328$source_name_ch1[x], split = " ")[[1]][3]
  if(status == "T") phen41328$source_name_ch1[x] <- "T" else phen41328$source_name_ch1[x] <- "N"
}

for(x in 1:nrow(phen74602)) {
  status <- strsplit(phen74602$source_name_ch1[x], split = " ")[[1]][1]
  if(status == "Normal") phen74602$source_name_ch1[x] <- "N" else phen74602$source_name_ch1[x] <- "T"
}

for(x in 1:nrow(phen142279)) {
  status <- phen142279$`tissue:ch1`[x]
  if(status == "adjacent normal") phen142279$`tissue:ch1`[x] <- "N" else phen142279$`tissue:ch1`[x] <- "T"
}

phen10950_ <- phen10950[,-2]
rownames(phen10950_) <- NULL
write.csv(phen10950_, "Results/mergeDatasets/phen10950.csv")

phen25070_ <- as.data.frame(cbind(sample=rownames(phen25070), group=phen25070[,1]))
write.csv(phen25070_, "Results/mergeDatasets/phen25070.csv")

phen41328_ <- as.data.frame(cbind(sample=rownames(phen41328), group=phen41328[,8]))
write.csv(phen41328_, "Results/mergeDatasets/phen41328.csv")

phen74602_ <- as.data.frame(cbind(sample=rownames(phen74602), group=phen74602[,8]))
write.csv(phen74602_, "Results/mergeDatasets/phen74602.csv")

phen142279_ <- as.data.frame(cbind(sample=rownames(phen142279), group=phen142279[,36]))
write.csv(phen142279_, "Results/mergeDatasets/phen142279.csv")
########################################


#### Load expression data ####
##############################
expr10950 <- as.data.frame(fread("Results/mergeDatasets/expr10950.csv"))
rownames(expr10950) <- expr10950$V1
expr10950 <- expr10950[,-1]

expr25070 <- as.data.frame(fread("Results/mergeDatasets/expr25070.csv"))
rownames(expr25070) <- expr25070$V1
expr25070 <- expr25070[,-1]

expr41328 <- as.data.frame(fread("Results/mergeDatasets/expr41328.csv"))
rownames(expr41328) <- expr41328$V1
expr41328 <- expr41328[,-1]

expr74602 <- as.data.frame(fread("Results/mergeDatasets/expr74602.csv"))
rownames(expr74602) <- expr74602$V1
expr74602 <- expr74602[,-1]

expr142279 <- as.data.frame(fread("Results/mergeDatasets/expr142279.csv"))
rownames(expr142279) <- expr142279$V1
expr142279 <- expr142279[,-1]

intersectSymbol <- Reduce(intersect, 
                          list(
                            rownames(expr10950),
                            rownames(expr25070),
                            rownames(expr41328),
                            rownames(expr74602),
                            rownames(expr142279)
                            )
                          )

expr10950 <- expr10950[intersectSymbol,]
expr25070 <- expr25070[intersectSymbol,]
expr41328 <- expr41328[intersectSymbol,]
expr74602 <- expr74602[intersectSymbol,]
expr142279 <- expr142279[intersectSymbol,]

expr10950 <- t(expr10950) %>% as.data.frame()
expr25070 <- t(expr25070) %>% as.data.frame()
expr41328 <- t(expr41328) %>% as.data.frame()
expr74602 <- t(expr74602) %>% as.data.frame()
expr142279 <- t(expr142279) %>% as.data.frame()
##############################

#### Load phenotypic data ####
##############################
phen10950 <- as.data.frame(fread("Results/mergeDatasets/phen10950.csv", header = TRUE))
phen10950 <- phen10950[,-1]
expr10950 <- cbind(sample=rownames(expr10950), expr10950)
expr10950 <- merge(phen10950, expr10950, by="sample")
rownames(expr10950) <- expr10950$sample
expr10950 <- expr10950[,-1]

phen25070 <- as.data.frame(fread("Results/mergeDatasets/phen25070.csv", header = TRUE))
phen25070 <- phen25070[,-1]
expr25070 <- cbind(sample=rownames(expr25070), expr25070)
expr25070 <- merge(phen25070, expr25070, by="sample")
rownames(expr25070) <- expr25070$sample
expr25070 <- expr25070[,-1]

phen41328 <- as.data.frame(fread("Results/mergeDatasets/phen41328.csv", header = TRUE))
phen41328 <- phen41328[,-1]
expr41328 <- cbind(sample=rownames(expr41328), expr41328)
expr41328 <- merge(phen41328, expr41328, by="sample")
rownames(expr41328) <- expr41328$sample
expr41328 <- expr41328[,-1]

phen74602 <- as.data.frame(fread("Results/mergeDatasets/phen74602.csv", header = TRUE))
phen74602 <- phen74602[,-1]
expr74602 <- cbind(sample=rownames(expr74602), expr74602)
expr74602 <- merge(phen74602, expr74602, by="sample")
rownames(expr74602) <- expr74602$sample
expr74602 <- expr74602[,-1]

phen142279 <- as.data.frame(fread("Results/mergeDatasets/phen142279.csv", header = TRUE))
phen142279 <- phen142279[,-1]
expr142279 <- cbind(sample=rownames(expr142279), expr142279)
expr142279 <- merge(phen142279, expr142279, by="sample")
rownames(expr142279) <- expr142279$sample
expr142279 <- expr142279[,-1]
##############################

#### Merging datasets ####
##########################
expr10950 <- cbind(batch=rep(1, nrow(expr10950)), expr10950)
expr25070 <- cbind(batch=rep(2, nrow(expr25070)), expr25070)
expr41328 <- cbind(batch=rep(3, nrow(expr41328)), expr41328)
expr74602 <- cbind(batch=rep(4, nrow(expr74602)), expr74602)
expr142279 <- cbind(batch=rep(5, nrow(expr142279)), expr142279)

mergedExpr <- rbind(expr10950, expr25070, expr41328, expr74602, expr142279)

mergedExpr$group <- as.factor(mergedExpr$group)
levels(mergedExpr$group) <- c(1,2)
mergedExpr$group <- as.numeric(mergedExpr$group)

for(x in 3:ncol(mergedExpr)) mergedExpr[,x] <- as.numeric(mergedExpr[,x])

write.csv(mergedExpr, "Results/mergeDatasets/mergeDatasets.csv")
##########################