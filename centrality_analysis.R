# load packages
library(igraph)
library(CINNA)
library(data.table)

# set working directory
setwd("C:/Users/ayoub/OneDrive/Desktop/")

#### Extracting giant component ####
giant.comp <- function(g) giant_component_extract(g)[[1]]

#### Creating PPI network ####
interactions <- read.delim("PPI/string_interactions_short.tsv", sep = "\t")

edges <- data.frame(node1=interactions$X.node1, node2=interactions$node2, score=interactions$combined_score)
nodes <- union(edges$node1, edges$node2)

ppi <- giant.comp(simplify(graph_from_data_frame(d = edges, vertices = nodes, directed = F)))

#### Centrality measure analysis ####
prop_cent <- proper_centralities(ppi)

# The following code takes about half hours to execute
system.time(calc_cent <- calculate_centralities(ppi, include = prop_cent[c(1:5,7:49)]))

pdf("PPI/centrality_analysis/pca_centralities.pdf")
pca_centralities(calc_cent)
dev.off()

closeness <- calculate_centralities(ppi, include = "Closeness centrality (Latora)")

ppiGenes <- sort(closeness$`Closeness centrality (Latora)`, decreasing = T)
ppiGenes <- rownames(as.data.frame(ppiGenes))

write.table(ppiGenes, file = "PPI/centrality_analysis/ppiGenes.txt", quote = F, row.names = F, col.names = F)
