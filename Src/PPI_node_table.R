library(data.table)
library(igraph)
library(CINNA)

setwd("D:/Sharif University/Master/Lessons/5. Fifth Term/Thesis/Res/")

interactions_node <- fread("PPI/string_interactions_short.tsv default node.csv")

#### Extracting giant component ####
giant.comp <- function(g) giant_component_extract(g)[[1]]

#### Creating PPI network ####
interactions <- fread("PPI/string_interactions_short.tsv")

edges <- data.frame(node1=interactions$`#node1`, 
                    node2=interactions$node2, 
                    score=interactions$combined_score)

nodes <- union(edges$node1, edges$node2)

ppi <- giant.comp(simplify(graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)))

measure <- "Diffusion Degree"
cnt <- calculate_centralities(ppi, include = measure)
cnt <- cnt$`Diffusion Degree`

cnt_measure_score <- data.frame(name=names(cnt), cnt_score=cnt)
interactions_node <- merge(interactions_node, cnt_measure_score, by="name")

upgenes <- read.table("limma/up.txt")
downgenes <- read.table("limma/down.txt")

colnames(upgenes) <- "gene"
colnames(downgenes) <- "gene"

up <- which(interactions_node$name %in% upgenes$gene)
down <- which(interactions_node$name %in% downgenes$gene)

interactions_node$regulated <- "down"
interactions_node$regulated[up] <- "up"
#interactions_node$regulated[down] <- "down"

write.table(interactions_node, file = "PPI/string_interactions_short.tsv default node.csv", sep = ",", row.names = F)
