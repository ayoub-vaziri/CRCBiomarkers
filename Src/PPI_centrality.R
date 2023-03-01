library(readr)
library(igraph)
library(CINNA)
library(ggraph)
library(graphlayouts)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/Res/")

#### Centrality measure analysis ####
#####################################

# Extracting giant component
giant.comp <- function(g) giant_component_extract(g)[[1]]

# Creating PPI network
interactions <- read_tsv("PPI/string_interactions_short.tsv")

edges <- data.frame(node1=interactions$`#node1`, 
                    node2=interactions$node2, 
                    score=interactions$combined_score)

nodes <- union(edges$node1, edges$node2)

ppi <- giant.comp(simplify(graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)))

# Centrality measure analysis
prop_cent <- proper_centralities(ppi)

calc_cent <- calculate_centralities(ppi, include = prop_cent[c(1:5,7:49)])

png("PPI/pca_centralities.png", height = 2200, width = 2600, res = 300)
options(repr.plot.width = 200, repr.plot.height = 100)
pca <- pca_centralities(calc_cent)
pca$labels$title <- ""
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
write.csv(pca_data, "PPI/pca_centrality_data.csv")

measure <- "Diffusion Degree"
cnt <- calculate_centralities(ppi, include = measure)
cnt <- sort(cnt$`Diffusion Degree`, decreasing = T)
cnt <- rownames(as.data.frame(cnt))

top20 <- head(cnt, 20)
write.table(top20, file = "PPI/top20.txt", quote = F, row.names = F, col.names = F)
#####################################