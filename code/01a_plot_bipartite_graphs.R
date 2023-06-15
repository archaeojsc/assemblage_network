require(ggraph)

# Plot bipartite assemblage network ---------------------------------------

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


# Plot random bipartite graph ---------------------------------------------

g_random_bpg_layout <- g_random_bpg %>% layout_as_bipartite()

g_random_bpg %>%
  ggraph(layout = g_random_bpg_layout) +
  geom_edge_link(edge_color = "darkgray", edge_alpha = 0.25) +
  geom_node_point(aes(color = type)) +
  ggtitle("Random Bipartite network graph")

