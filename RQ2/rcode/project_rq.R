library(dplyr)
library(tidyr)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggplot2)
# optional for enrichment (uses web API)
install_if_missing <- function(pk){ if(!requireNamespace(pk, quietly=TRUE)) install.packages(pk) }
install_if_missing("gprofiler2")
library(gprofiler2)

load("results/final_networks.RData")  # bipartite_graph, ppi_graph, combined_graph, drug_target_oncology, oncology_map

# nodes/edges counts
cat("Bipartite nodes:", vcount(bipartite_graph), "edges:", ecount(bipartite_graph), "\n")
cat("PPI nodes:", vcount(ppi_graph), "edges:", ecount(ppi_graph), "\n")
cat("Combined nodes:", vcount(combined_graph), "edges:", ecount(combined_graph), "\n")

# helper: convert igraph to tidygraph
bg_tbl <- as_tbl_graph(bipartite_graph)
ppi_tbl <- as_tbl_graph(ppi_graph)
comb_tbl <- as_tbl_graph(combined_graph)

# RQ1

# create drug-protein edge table (drug_target_oncology is from ChEMBL; if you used full dataset it's already big)
dt_edges <- drug_target_oncology %>%
  mutate(Drug = toupper(DrugName), Protein = toupper(UniProt)) %>%
  select(Drug, Protein) %>%
  distinct()

# build bipartite incidence matrix (sparse)
library(Matrix)
drugs <- sort(unique(dt_edges$Drug))
prots <- sort(unique(dt_edges$Protein))
i <- match(dt_edges$Drug, drugs)
j <- match(dt_edges$Protein, prots)
M <- sparseMatrix(i = i, j = j, dims = c(length(drugs), length(prots)))
# project to drug-drug co-target matrix (count of shared proteins)
drug_drug_mat <- M %*% t(M)
diag(drug_drug_mat) <- 0

# create igraph object for drug projection (only keep edges with weight >= 1 or set threshold)
# convert sparse matrix to edge list
library(Matrix)
drug_pairs <- which(drug_drug_mat >= 1, arr.ind = TRUE)
weights <- drug_drug_mat[drug_pairs]
edges_df <- tibble(
  from = drugs[drug_pairs[,1]],
  to   = drugs[drug_pairs[,2]],
  weight = as.numeric(weights)
) %>% filter(as.character(from) < as.character(to))  # keep unique undirected

drug_proj_graph <- graph_from_data_frame(edges_df, directed = FALSE, vertices = tibble(name = drugs))
# community detection (Louvain)
drug_comm <- cluster_louvain(drug_proj_graph, weights = E(drug_proj_graph)$weight)
V(drug_proj_graph)$comm <- membership(drug_comm)
table(V(drug_proj_graph)$comm) %>% head()
# save results
drug_comm_df <- tibble(Drug = V(drug_proj_graph)$name, DrugCommunity = V(drug_proj_graph)$comm)
write.csv(drug_comm_df, "results/drug_communities_projection.csv", row.names = FALSE)


# VISUALIZE

# Visualize the drug projection (keep top drugs by degree or community)
deg_drug <- degree(drug_proj_graph)
top_drugs <- names(sort(deg_drug, decreasing = TRUE))[1:80]  # pick top 80 hubs for figure
subg <- induced_subgraph(drug_proj_graph, vids = which(V(drug_proj_graph)$name %in% top_drugs))

library(ggraph)
set.seed(42)
ggraph(subg, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.3) +
  geom_node_point(aes(color = as.factor(comm)), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  labs(color = "Drug comm")
ggsave("results/plot_drug_projection_top80.png", width = 14, height = 10, bg = 'white')


# RQ2

# create mapping Drug -> DrugCommunity (from Approach A)
drug_comm_df <- read.csv("results/drug_communities_projection.csv", stringsAsFactors = FALSE)

# build drug->protein list from drug_target_oncology; ensure Drug names match drug_comm_df
dt <- drug_target_oncology %>% mutate(Drug = toupper(DrugName), Protein = toupper(UniProt)) %>% select(Drug, Protein) %>% distinct()
dt <- dt %>% left_join(drug_comm_df, by = c("Drug" = "Drug"))

# for each protein, count distinct drug communities that target it
protein_comm_stats <- dt %>%
  group_by(Protein) %>%
  summarize(
    n_drugs = n_distinct(Drug),
    n_drug_comms = n_distinct(DrugCommunity, na.rm = TRUE),
    drug_comms = paste0(sort(unique(DrugCommunity[!is.na(DrugCommunity)])), collapse = ";")
  ) %>% arrange(desc(n_drug_comms), desc(n_drugs))

# compute betweenness centrality on PPI for proteins (ensure node IDs match: PPI graph must have Protein IDs)
# If PPI nodes are gene symbols and your protein list is UniProt, you need mapping. We'll attempt betweenness on combined_graph as fallback:
bet_protein <- betweenness(ppi_graph, normalized = TRUE)  # by gene symbol nodes
bet_df <- tibble(node = names(bet_protein), betweenness = as.numeric(bet_protein))

# Try to join protein_comm_stats with bet_df by matching Protein->node names.
bridge_candidates <- protein_comm_stats %>%
  left_join(bet_df, by = c("Protein" = "node")) %>%
  arrange(desc(n_drug_comms), desc(betweenness, na.rm = TRUE)) %>%
  relocate(Protein, n_drugs, n_drug_comms, betweenness)

write.csv(bridge_candidates, "results/protein_bridge_candidates.csv", row.names = FALSE)
head(bridge_candidates, 30)

# Example: enrichment for top protein community (get proteins in a PPI community)
# pick a community index from ppi_comm (if computed)
if (exists("ppi_comm")) {
  top_pcomm <- which.max(sizes(ppi_comm))  # largest community as example
  p_nodes <- names(which(membership(ppi_comm) == top_pcomm))
  # run gprofiler (specify organism = "hsapiens")
  gost_res <- gost(p_nodes, organism = "hsapiens", sources = c("GO:BP","KEGG","REAC"))
  gost_res$result %>% head()
  # save results
  write.csv(gost_res$result, paste0("results/enrichment_ppi_comm_", top_pcomm, ".csv"), row.names = FALSE)
}

# Enrichment for top bridge proteins
top_bridges <- bridge_candidates %>% filter(!is.na(betweenness)) %>% arrange(desc(n_drug_comms), desc(betweenness)) %>% head(50)
gost_bridges <- gost(top_bridges$Protein, organism = "hsapiens", sources = c("GO:BP","KEGG","REAC"))
write.csv(gost_bridges$result, "results/enrichment_top_bridge_proteins.csv", row.names = FALSE)

# VISUALIZE

# PPI subgraph: pick a community or top proteins
if (exists("ppi_comm")) {
  comm_id <- which.max(sizes(ppi_comm))
  p_nodes <- names(which(membership(ppi_comm) == comm_id))
  subppi <- induced_subgraph(ppi_graph, vids = which(V(ppi_graph)$name %in% p_nodes))
  # overlay target annotation
  V(subppi)$is_target <- V(subppi)$name %in% dt$Protein
  # compute betweenness for subgraph
  bsub <- betweenness(subppi, normalized = TRUE)
  V(subppi)$bet <- bsub
  # plot
  ggraph(subppi, layout = "fr") +
    geom_edge_link(alpha = 0.2) +
    geom_node_point(aes(size = bet, color = is_target)) +
    geom_node_text(aes(label = ifelse(is_target, name, "")), repel = TRUE, size = 3) +
    theme_void()
  ggsave("results/ppi_community_snapshot.png", width=10, height=8)
}

bridge_table <- bridge_candidates %>%
  left_join(
    dt %>% group_by(Protein) %>% summarize(drugs = paste(unique(Drug), collapse = "; "), n_drugs = n_distinct(Drug)),
    by = "Protein"
  ) %>%
  arrange(desc(n_drug_comms), desc(betweenness)) %>%
  select(Protein, n_drug_comms, n_drugs, betweenness, drug_comms, drugs) %>%
  distinct()
write.csv(bridge_table, "results/top_bridge_table.csv", row.names = FALSE)

gc <- decompose.graph(ppi_graph)[[1]]