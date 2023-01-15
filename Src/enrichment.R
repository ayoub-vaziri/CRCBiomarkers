library(data.table)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/")

symbol <- fread("Res/PPI/DEGsRGs.txt", header = F)$V1

#### GO enrichment analysis ####
################################
ego <- enrichGO(gene = symbol, 
               OrgDb = "org.Hs.eg.db", 
               ont = "ALL",
               keyType = "SYMBOL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)

png("Res/enrichment/GOenrichment.png", height = 2200, width = 2200, res = 300)
barplot(ego, 
        split = "ONTOLOGY", 
        showCategory = 10, 
        font.size = 8.5, 
        color = "qvalue", 
        label_format = 60) +
  facet_grid(ONTOLOGY~., scale = "free")
dev.off()
################################

entrezid <- mapIds(x = org.Hs.eg.db, keys = symbol, column = "ENTREZID", keytype = "SYMBOL") %>% 
  as.data.frame()

#### KEGG enrichment analysis ####
##################################
ekegg <- enrichKEGG(gene = entrezid$.,
           organism = "hsa",
           keyType = "kegg",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           qvalueCutoff = 0.05)

png("Res/enrichment/KEGGenrichment.png", height = 2000, width = 2000, res = 300)
dotplot(ekegg, 
        showCategory = 10, 
        font.size = 8, 
        color = "qvalue", 
        label_format = 70)
dev.off()
##################################

#### DO enrichment analysis ####
################################
edo <- enrichDO(gene = entrezid$.,
                ont = "DO",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05
)

png("Res/enrichment/DOenrichment.png", height = 2000, width = 2800, res = 300)
barplot(edo, 
        showCategory = 10, 
        font.size = 10, 
        color = "qvalue", 
        label_format = 70,
        )
dev.off()
################################

#### Gene-concept network for KEGG ####
#######################################
ekeggx <- setReadable(x = ekegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetkegg <- cnetplot(ekeggx, 
                     circular = TRUE, 
                     colorEdge = TRUE,
                     showCategory = 10,
                     cex_label_category = 0.72,
                     cex_label_gene = 0.6
                     )

png("Res/enrichment/KEGGcnetplot.png", height = 2000, width = 3200, res = 300)
cnetkegg
dev.off()
#######################################

#### Gene-concept network for GO ####
#####################################
egox <- setReadable(x = ego, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetgo <- cnetplot(egox, 
                   circular = TRUE, 
                   colorEdge = TRUE, 
                   showCategory = 10,
                   cex_label_category = 0.72,
                   cex_label_gene = 0.6)

png("Res/enrichment/GOcnetplot.png", height = 2000, width = 3400, res = 300)
cnetgo
dev.off()
#####################################

#### DO enrichment analysis ####
################################
edox <- setReadable(x = edo, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetdo <- cnetplot(edox, 
                   circular = TRUE, 
                   colorEdge = TRUE, 
                   showCategory = 5,
                   cex_label_category = 0.8,
                   cex_label_gene = 0.6)

png("Res/enrichment/DOcnetplot.png", height = 2000, width = 3000, res = 300)
cnetdo
dev.off()
################################