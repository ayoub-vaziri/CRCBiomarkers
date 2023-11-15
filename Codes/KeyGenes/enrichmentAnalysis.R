library(data.table)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

symbol <- fread("Res/differentialExpressionAnalysis/updown.txt", header = F)$V1

entrezid <- mapIds(x = org.Hs.eg.db, keys = symbol, column = "ENTREZID", keytype = "SYMBOL") %>% 
  as.data.frame() %>% na.omit()

colnames(entrezid) <- "ENTREZID"

#### GO enrichment analysis ####
################################
cc <- enrichGO(gene = symbol, 
               OrgDb = "org.Hs.eg.db", 
               ont = "CC",
               keyType = "SYMBOL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)

mf <- enrichGO(gene = symbol, 
               OrgDb = "org.Hs.eg.db", 
               ont = "MF",
               keyType = "SYMBOL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)

bp <- enrichGO(gene = symbol, 
               OrgDb = "org.Hs.eg.db", 
               ont = "BP",
               keyType = "SYMBOL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)

png("Res/enrichmentAnalysis/CC.png", height = 2000, width = 2600, res = 300)
dotplot(cc, 
        # split = "ONTOLOGY", 
        showCategory = 10, 
        font.size = 14, 
        color = "qvalue", 
        label_format = 60) 
  # facet_grid(cc~., scale = "free")
dev.off()

png("Res/enrichmentAnalysis/MF.png", height = 2000, width = 2600, res = 300)
dotplot(mf, 
        # split = "ONTOLOGY", 
        showCategory = 10, 
        font.size = 14, 
        color = "qvalue", 
        label_format = 60) 
  # facet_grid(ONTOLOGY~., scale = "free")
dev.off()

png("Res/enrichmentAnalysis/BP.png", height = 2000, width = 2600, res = 300)
dotplot(bp, 
        # split = "ONTOLOGY", 
        showCategory = 10, 
        font.size = 14, 
        color = "qvalue", 
        label_format = 60) 
  # facet_grid(ONTOLOGY~., scale = "free")
dev.off()
################################

#### KEGG enrichment analysis ####
##################################
ekegg <- enrichKEGG(gene = entrezid$ENTREZID,
           organism = "hsa",
           keyType = "kegg",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH")

png("Res/enrichmentAnalysis/KEGG.png", height = 2000, width = 2600, res = 300)
dotplot(ekegg, 
        showCategory = 10, 
        font.size = 14, 
        color = "qvalue"
        )
dev.off()
##################################