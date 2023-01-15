library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(annotate)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

#### Load data ####
###################
gse25070 <- getGEO("GSE25070", destdir = "Data/")
gse25071 <- getGEO("GSE25071", destdir = "Data/")
gse44076 <- getGEO("GSE44076", destdir = "Data/")
gse106582 <- getGEO("GSE106582", destdir = "Data/")
gse110223 <- getGEO("GSE110223", destdir = "Data/")
gse110224 <- getGEO("GSE110224", destdir = "Data/")

gse25070 <- gse25070[[1]]
gse25071 <- gse25071[[1]]
gse44076 <- gse44076[[1]]
gse106582 <- gse106582[[1]]
gse110223 <- gse110223[[1]]
gse110224 <- gse110224[[1]]

expr25070 <- exprs(gse25070)
expr25071 <- exprs(gse25071)
expr44076 <- exprs(gse44076)
expr106582 <- exprs(gse106582)
expr110223 <- exprs(gse110223)
expr110224 <- exprs(gse110224)

phen25070 <- pData(gse25070)
phen25071 <- pData(gse25071)
phen44076 <- pData(gse44076)
phen106582 <- pData(gse106582)
phen110223 <- pData(gse110223)
phen110224 <- pData(gse110224)

feat25070 <- fData(gse25070)
feat25071 <- fData(gse25071)
feat44076 <- fData(gse44076)
feat106582 <- fData(gse106582)
feat110223 <- fData(gse110223)
feat110224 <- fData(gse110224)
###################

#### Process data ####
######################
entrez25070 <- as.character(unique(feat25070$Entrez_Gene_ID))
entrez25071 <- unique(unlist(strsplit(feat25071$GENE, split = ',')))
entrez44076 <- na.omit(unique(unlist(strsplit(feat44076$`Entrez Gene`, split = " /// "))))
entrez106582 <- as.character(unique(feat106582$Entrez_Gene_ID))
entrez110223 <- unique(unlist(strsplit(feat110223$ENTREZ_GENE_ID, split = " /// ")))
entrez110224 <- unique(unlist(strsplit(feat110224$ENTREZ_GENE_ID, split = " /// ")))

intersectEntrez <- Reduce(intersect, list(entrez25070, entrez25071, entrez44076, entrez106582, entrez110223, entrez110224))

getSymbol <- function(entrez_id) {
  entrez_id <- as.character(entrez_id)
  return(lookUp(entrez_id, 'org.Hs.eg.db', 'SYMBOL')[[1]])
}

entrez.to.symbol <- function(feat, spliter = NULL, symbolColName, entrez_col = NULL) {
  
  if(is.null(entrez_col))
    entrez_col <- grep("Entrez.*", colnames(feat), ignore.case = TRUE)
  feat <- cbind(feat, Entrez_ID=feat[,entrez_col])
  feat <- cbind(feat, symbol=feat[,symbolColName])
  
  for(eid in intersectEntrez) {
    rowx <- grep(eid, feat$Entrez_ID)
    ids <- lapply(rowx, function(x) {
      if(is.null(spliter)) 
        return(nchar(as.character(feat$Entrez_ID[x])) == nchar(eid))
      else {
        ss <- unlist(strsplit(as.character(feat$Entrez_ID[x]), split = spliter))
        stat <- which(nchar(ss) == nchar(eid) & ss == eid)
        if(identical(stat, integer(0))) return(FALSE) else return(TRUE) 
      }
    }) %>% unlist()

    if(!is.null(ids) & !identical(rowx[ids[1]], integer(0))) {
      ids <- which(ids)
      feat$Entrez_ID[rowx[ids[1]]] <- eid
    }
  }
  
  for(x in 1:nrow(feat)) {
    if(is.null(spliter))
      entrez_id <- as.character(feat$Entrez_ID[x])
    else
      entrez_id <- strsplit(as.character(feat$Entrez_ID[x]), split = spliter)[[1]][1]
    if(!is.na(entrez_id) & !is.null(entrez_id) & !identical(entrez_id, integer(0)) & !identical(entrez_id, "")) {
      gene_symbol <- getSymbol(entrez_id)
      if(!is.na(gene_symbol) & !is.null(gene_symbol) & !identical(gene_symbol, integer(0)))
        feat$symbol[x] <- gene_symbol
    }  
  }
  
  if(!is.null(spliter)) {
    for(x in 1:nrow(feat)) {
      ss <- strsplit(as.character(feat$symbol[x]), split = spliter)[[1]]
      if(length(ss) > 1)
        feat$symbol[x] <- ss[1]
    }
  }
  
  return(feat)
}

feat25070 <- entrez.to.symbol(feat25070, symbolColName = "Symbol")
feat25071 <- entrez.to.symbol(feat25071, symbolColName = "Gene Symbol", spliter = ',', entrez_col = "GENE")
feat44076 <- entrez.to.symbol(feat44076, spliter = " /// ", symbolColName = "Gene Symbol")
feat106582 <- entrez.to.symbol(feat106582, symbolColName = "Symbol")
feat110223 <- entrez.to.symbol(feat110223, spliter = " /// ", symbolColName = "Gene Symbol")
feat110224 <- entrez.to.symbol(feat110224, spliter = " /// ", symbolColName = "Gene Symbol")

expr25070 <- cbind(ID=rownames(expr25070), expr25070) %>% as.data.frame()
expr25071 <- cbind(ID=rownames(expr25071), expr25071) %>% as.data.frame()
expr44076 <- cbind(ID=rownames(expr44076), expr44076) %>% as.data.frame()
expr106582 <- cbind(ID=rownames(expr106582), expr106582) %>% as.data.frame()
expr110223 <- cbind(ID=rownames(expr110223), expr110223) %>% as.data.frame()
expr110224 <- cbind(ID=rownames(expr110224), expr110224) %>% as.data.frame()

feat25070 <- feat25070[, c("ID", "symbol")]
feat25071 <- feat25071[, c("ID", "symbol")]
feat44076 <- feat44076[, c("ID", "symbol")]
feat106582 <- feat106582[, c("ID", "symbol")]
feat110223 <- feat110223[, c("ID", "symbol")]
feat110224 <- feat110224[, c("ID", "symbol")]

expr25070 <- merge(feat25070, expr25070, by = "ID")
expr25071 <- merge(feat25071, expr25071, by = "ID")
expr44076 <- merge(feat44076, expr44076, by = "ID")
expr106582 <- merge(feat106582, expr106582, by = "ID")
expr110223 <- merge(feat110223, expr110223, by = "ID")
expr110224 <- merge(feat110224, expr110224, by = "ID")

expr25070 <- expr25070[which(expr25070$symbol != ""),]
expr25071 <- expr25071[which(expr25071$symbol != ""),]
expr44076 <- expr44076[which(expr44076$symbol != ""),]
expr106582 <- expr106582[which(expr106582$symbol != ""),]
expr110223 <- expr110223[which(expr110223$symbol != ""),]
expr110224 <- expr110224[which(expr110224$symbol != ""),]

expr25070 <- expr25070[!duplicated(expr25070$symbol),]
expr25071 <- expr25071[!duplicated(expr25071$symbol),]
expr44076 <- expr44076[!duplicated(expr44076$symbol),]
expr106582 <- expr106582[!duplicated(expr106582$symbol),]
expr110223 <- expr110223[!duplicated(expr110223$symbol),]
expr110224 <- expr110224[!duplicated(expr110224$symbol),]

rownames(expr25070) <- expr25070$symbol
rownames(expr25071) <- expr25071$symbol
rownames(expr44076) <- expr44076$symbol
rownames(expr106582) <- expr106582$symbol
rownames(expr110223) <- expr110223$symbol
rownames(expr110224) <- expr110224$symbol

expr25070 <- expr25070[,-c(1,2)]
expr25071 <- expr25071[, -c(1,2)]
expr44076 <- expr44076[,-c(1,2)]
expr106582 <- expr106582[,-c(1,2)]
expr110223 <- expr110223[,-c(1,2)]
expr110224 <- expr110224[,-c(1,2)]

expr25070 <- na.omit(expr25070)
expr25071 <- na.omit(expr25071)
expr44076 <- na.omit(expr44076)
expr106582 <- na.omit(expr106582)
expr110223 <- na.omit(expr110223)
expr110224 <- na.omit(expr110224)

write.csv(expr25070, "Res/exprBatchBioData/expr25070.csv")
write.csv(expr25071, "Res/exprBatchBioData/expr25071.csv")
write.csv(expr44076, "Res/exprBatchBioData/expr44076.csv")
write.csv(expr106582, "Res/exprBatchBioData/expr106582.csv")
write.csv(expr110223, "Res/exprBatchBioData/expr110223.csv")
write.csv(expr110224, "Res/exprBatchBioData/expr110224.csv")

expr25070 <- as.data.frame(fread("Res/exprBatchBioData/expr25070.csv"))
rownames(expr25070) <- expr25070$V1
expr25070 <- expr25070[,-1]

expr25071 <- as.data.frame(fread("Res/exprBatchBioData/expr25071.csv"))
rownames(expr25071) <- expr25071$V1
expr25071 <- expr25071[,-1]

expr44076 <- as.data.frame(fread("Res/exprBatchBioData/expr44076.csv"))
rownames(expr44076) <- expr44076$V1
expr44076 <- expr44076[,-1]

expr106582 <- as.data.frame(fread("Res/exprBatchBioData/expr106582.csv"))
rownames(expr106582) <- expr106582$V1
expr106582 <- expr106582[,-1]

expr110223 <- as.data.frame(fread("Res/exprBatchBioData/expr110223.csv"))
rownames(expr110223) <- expr110223$V1
expr110223 <- expr110223[,-1]

expr110224 <- as.data.frame(fread("Res/exprBatchBioData/expr110224.csv"))
rownames(expr110224) <- expr110224$V1
expr110224 <- expr110224[,-1]

intersectSymbol <- Reduce(intersect, list(rownames(expr25070), rownames(expr25071), rownames(expr44076), rownames(expr106582), rownames(expr110223), rownames(expr110224)))

expr25070 <- expr25070[intersectSymbol,]
expr25071 <- expr25071[intersectSymbol,]
expr44076 <- expr44076[intersectSymbol,]
expr106582 <- expr106582[intersectSymbol,]
expr110223 <- expr110223[intersectSymbol,]
expr110224 <- expr110224[intersectSymbol,]

expr25070 <- t(expr25070) %>% as.data.frame()
expr25071 <- t(expr25071) %>% as.data.frame()
expr44076 <- t(expr44076) %>% as.data.frame()
expr106582 <- t(expr106582) %>% as.data.frame()
expr110223 <- t(expr110223) %>% as.data.frame()
expr110224 <- t(expr110224) %>% as.data.frame()

exprMerged <- rbind(expr25070, expr25071, expr44076, expr106582, expr110223, expr110224)

for(x in 1:ncol(exprMerged)) exprMerged[,x] <- as.numeric(exprMerged[,x])

write.csv(exprMerged, "Res/exprBatchBioData/exprColorectal.csv")
######################

#### Identifying batches ####
#############################
batch25070 <- rep(1, nrow(expr25070))
batch25071 <- rep(2, nrow(expr25071))
batch44076 <- rep(3, nrow(expr44076))
batch106582 <- rep(4, nrow(expr106582))
batch110223 <- rep(5, nrow(expr110223))
batch110224 <- rep(6, nrow(expr110224))

batches <- c(batch25070, batch25071, batch44076, batch106582, batch110223, batch110224)

write.table(batches, file = "Res/exprBatchBioData/batchColorectal.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
#############################

#### Identifying biological classes ####
########################################
for(x in 1:nrow(phen25070)) {
  status <- strsplit(phen25070$title[x], split = " ")[[1]][1]
  if(status == "Colorectal") phen25070$title[x] <- "T" else phen25070$title[x] <- "N"
}

for(x in 1:nrow(phen25071)) {
  status <- strsplit(phen25071$source_name_ch1[x], split = " ")[[1]][1]
  if(status == "colorectal") phen25071$source_name_ch1[x] <- "T" else phen25071$source_name_ch1[x] <- "N"
}

for(x in 1:nrow(phen44076)) {
  status <- strsplit(phen44076$title[x], split = " ")[[1]][1]
  if(status == "Tumor") phen44076$title[x] <- "T" else phen44076$title[x] <- "N"
}

for(x in 1:nrow(phen106582)) {
  status <- strsplit(phen106582$source_name_ch1[x], split = " ")[[1]][4]
  if(status == "tumor") phen106582$source_name_ch1[x] <- "T" else phen106582$source_name_ch1[x] <- "N"
}

for(x in 1:nrow(phen110223)) {
  status <- strsplit(phen110223$title[x], split = "_")[[1]][3]
  if(status == "cancer") phen110223$title[x] <- "T" else phen110223$title[x] <- "N"
}

for(x in 1:nrow(phen110224)) {
  status <- strsplit(phen110224$title[x], split = "_")[[1]][3]
  if(substr(status,1,6) == "cancer") phen110224$title[x] <- "T" else phen110224$title[x] <- "N"
}

bio <- c(phen25070$title, phen25071$source_name_ch1, phen44076$title, phen106582$source_name_ch1, phen110223$title, phen110224$title)
bio <- as.factor(bio)
levels(bio) <- c(1,2)
bio <- as.numeric(bio)

write.table(bio, file = "Res/exprBatchBioData/bioColorectal.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
########################################