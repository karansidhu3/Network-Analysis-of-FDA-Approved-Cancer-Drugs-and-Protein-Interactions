library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

# Convert bipartite_graph to tbl_graph for tidy plotting (optional)
bg_tbl <- as_tbl_graph(bipartite_graph)

# Assign type labels for clarity
V(bg_tbl)$label <- V(bg_tbl)$name
V(bg_tbl)$type_label <- ifelse(V(bg_tbl)$type == 0, "Drug", "Protein")

# Optional: detect communities on the drug projection for color coding
# Create drug-protein edges
drug_edges <- bipartite_graph %>% 
  bipartite_projection(which = "true")  # true = first mode = drugs
drug_comm <- cluster_louvain(drug_edges)
V(bipartite_graph)$comm <- NA
V(bipartite_graph)$comm[V(bipartite_graph)$type == 0] <- membership(drug_comm)

# Plot bipartite graph
set.seed(42)
p <- ggraph(bg_tbl, layout = "bipartite") +
  geom_edge_link(alpha = 0.3, color = "grey70") +
  geom_node_point(aes(color = factor(comm), shape = type_label), size = 5) +
  geom_node_text(aes(label = ifelse(type_label == "Drug", label, "")), 
                 repel = TRUE, size = 3, max.overlaps = 50) +
  scale_shape_manual(values = c(16, 17)) +  # circles for drugs, triangles for proteins
  scale_color_viridis_d(option = "C", na.value = "grey50") +
  theme_void() +
  theme(legend.position = "right") +
  labs(color = "Drug Community", shape = "Node Type")

# Save as PNG
ggsave("results/bipartite_graph_clean.png", p, width = 14, height = 10, bg = "white")
