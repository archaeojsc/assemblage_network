require(tidyverse)
require(igraph)
require(ggraph)


# Create un-weighted bipartite graph --------------------------------------

g_assemblages_bpg <-
  graph_from_data_frame(unique.data.frame(select(dat, LEVEL_ID, CODE)),
                        directed = FALSE)


# Assign vertex types -----------------------------------------------------

V(g_assemblages_bpg)$type <-
  bipartite_mapping(g_assemblages_bpg)$type


# Create incidence matrix from bipartite graph -----------------------------

g_assemblages_bpg_inc <- as_incidence_matrix(g_assemblages_bpg)


# Plot bipartite graph ----------------------------------------------------

g_assemblages_bpg_layout <-
  g_assemblages_bpg %>% layout_as_bipartite()

g_assemblages_bpg %>%
  ggraph(layout = g_assemblages_bpg_layout) +
  geom_edge_link(edge_color = "darkgray", edge_alpha = 0.25) +
  geom_node_point(aes(color = type)) +
  scale_color_manual(
    values = c("darkgreen", "darkblue"),
    name = "Node Type",
    breaks = c(FALSE, TRUE),
    labels = c("Provenience", "Artifact")
  ) +
  ggtitle("Bipartite network of Provenience and Artifact Type")


# Create random bipartite graph (for comparisons) -------------------------

g_random_bpg <-
  sample_bipartite(
    n1 = 152,
    n2 = 233,
    type = "gnm",
    m = 2240,
    directed = FALSE
  )

g_random_bpg_inc <- as_incidence_matrix(g_random_bpg)

g_random_bpg_layout <- g_random_bpg %>% layout_as_bipartite()

g_random_bpg %>%
  ggraph(layout = g_random_bpg_layout) +
  geom_edge_link(edge_color = "darkgray", edge_alpha = 0.25) +
  geom_node_point(aes(color = type)) +
  ggtitle("Random Bipartite network graph")

