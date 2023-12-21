# g_prov_projection %>%
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "darkgreen") +
#   ggtitle("Network of Proveniences")
# 
# g_artifact_projection %>%
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "darkblue") +
#   ggtitle("Network of Artifact Types")

prov_strength <-
  data.frame(ssoc = strength(g_prov_oc),
             jacc = strength(g_prov_jacc),
             sd = strength(g_prov_sd))

prov_strength %>% stack() %>%
  ggplot(aes(x = values)) +
  geom_density(fill = "darkgreen",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Node Strength for Proveniences Graph")

rand_t_strength <-
  data.frame(ssoc = strength(g_rand_t_oc),
             jacc = strength(g_rand_t_jacc),
             sd = strength(g_rand_t_sd))

rand_t_strength %>% stack() %>%
  ggplot(aes(x = values)) +
  geom_density(fill = "darkgreen",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Node Strength for Random 'Provenience' Graph")

artifact_strength <-
  data.frame(ssoc = strength(g_artifact_oc),
             jacc = strength(g_artifact_jacc),
             sd = strength(g_artifact_sd))

artifact_strength %>% stack() %>%
  ggplot(aes(x = values)) +
  geom_density(fill = "darkblue",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Node Strength for Artifacts Graph")

rand_b_strength <-
  data.frame(ssoc = strength(g_rand_b_oc),
             jacc = strength(g_rand_b_jacc),
             sd = strength(g_rand_b_sd))

rand_b_strength %>% stack() %>%
  ggplot(aes(x = values)) +
  geom_density(fill = "darkblue",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Node Strength for Random 'Artifacts' Graph")


rand_t_degree <-
  data.frame(ssoc = degree(g_rand_t_oc),
             jacc = degree(g_rand_t_jacc),
             sd = degree(g_rand_t_sd))

rand_t_degree %>% stack() %>%
  ggplot(aes(x = values)) +
  geom_density(fill = "darkgreen",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Node Degree for Random 'Provenience' Graph")

rand_t_edges <-
  data.frame(ssoc = E(g_rand_t_oc)$weight,
             jacc = E(g_rand_t_jacc)$weight,
             sd = E(g_rand_t_sd)$weight)

rand_t_edges %>% stack() %>%
  ggplot(aes(x = values)) +
  geom_density(fill = "darkgreen",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Edge Weights for Random 'Provenience' Graph")

