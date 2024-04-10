library(data.table)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Set the current working directory to the project path
setwd(project_path)

symbol <- fread("Results/differentialExpressionAnalysis/updown.txt", header = F)$V1

entrezid <- mapIds(x = org.Hs.eg.db, keys = symbol, column = "ENTREZID", keytype = "SYMBOL") %>% 
  as.data.frame() %>% na.omit()

colnames(entrezid) <- "ENTREZID"

#### GO enrichment analysis ####
################################
all <- enrichGO(gene = symbol, 
                OrgDb = "org.Hs.eg.db", 
                ont = "ALL",
                keyType = "SYMBOL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)

png("Results/enrichmentAnalysis/ALL.png", height = 2000, width = 2600, res = 300)
barplot(all, 
        split = "ONTOLOGY", 
        showCategory = 10, 
        font.size = 10, 
        color = "qvalue", 
        label_format = 60) + 
  facet_grid(ONTOLOGY~., scale = "free")
dev.off()
################################

#### KEGG enrichment analysis ####
##################################
ekegg <- enrichKEGG(gene = entrezid$ENTREZID,
           organism = "hsa",
           keyType = "kegg",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH")

png("Results/enrichmentAnalysis/KEGG.png", height = 1700, width = 2000, res = 300)
dotplot(ekegg, 
        showCategory = 10, 
        font.size = 12, 
        color = "qvalue"
        )
dev.off()
##################################