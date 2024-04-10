library(readr)
library(igraph)
library(CINNA)
library(ggraph)
library(graphlayouts)

# Set the current working directory to the project path
setwd(project_path)

#### Centrality measure analysis ####
#####################################

# Extracting giant component
giant.comp <- function(g) giant_component_extract(g)[[1]]

# Creating PPI network
interactions <- read_tsv("Results/centralityAnalysis/string_interactions_short.tsv")

edges <- data.frame(node1=interactions$`#node1`, 
                    node2=interactions$node2, 
                    score=interactions$combined_score)

nodes <- union(edges$node1, edges$node2)

write.csv(edges, "Results/centralityAnalysis//edges.csv")
write.csv(as.data.frame(nodes), "Results/centralityAnalysis//nodes.csv")

ppi <- giant.comp(simplify(graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)))

# Centrality measure analysis
prop_cent <- proper_centralities(ppi)

calc_cent <- calculate_centralities(ppi, include = prop_cent[c(1:5, 7:17, 19:26, 28, 32, 34:36, 40:42, 44:49)])

png("Results/centralityAnalysis//pca_centralities.png", height = 2200, width = 2600, res = 300)
options(repr.plot.width = 200, repr.plot.height = 100)
pca <- pca_centralities(calc_cent)
pca + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 9),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.margin = margin(rep(15, 4)))
dev.off()

pca_data <- pca$data
pca_data <- pca_data[order(pca_data$contrib, decreasing = TRUE),]
rownames(pca_data) <- NULL
colnames(pca_data) <- c("Centrality measure", "Contribution")
write.csv(pca_data, "Results/centralityAnalysis/pca_centrality_data.csv")

measure <- "Closeness centrality (Latora)"
cnt <- calculate_centralities(ppi, include = measure)
cnt <- cnt$`Closeness centrality (Latora)`
cnt <- cnt[cnt >= mean(cnt)]
cnt <- rownames(as.data.frame(cnt))

write.table(cnt, file = "Results/centralityAnalysis/essentialNodes.txt", quote = F, row.names = F, col.names = F)
#####################################