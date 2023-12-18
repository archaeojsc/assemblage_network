


# Run previous scripts ----------------------------------------------------

source("code/03_similarity_and_adjacency.R", echo = TRUE)

# Project with igraph  ----------------------------------------------------


# Weighted projection for each node type

g_assemblages_proj <-
  bipartite_projection(g_assemblages_bpg, multiplicity = TRUE)

g_prov_projection <- g_assemblages_proj$proj1

g_artifact_projection <- g_assemblages_proj$proj2

g_random_proj <-
  bipartite.projection(g_random_bpg, multiplicity = TRUE)

g_rand_t_projection <- g_random_proj$proj1

g_rand_b_projection <- g_random_proj$proj2


# Project sub-graphs from overlap coefficient ----------------------------


g_prov_oc <-
  graph_from_adjacency_matrix(sim_prov_oc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)


g_artifact_oc <-
  graph_from_adjacency_matrix(sim_artifact_oc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)


g_rand_t_oc <-
  graph_from_adjacency_matrix(sim_rand_t_oc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

g_rand_b_oc <-
  graph_from_adjacency_matrix(sim_rand_b_oc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)


# Project with Sorenson-Dice ---------------------------------------------

g_prov_sd <-
  graph_from_adjacency_matrix(sim_prov_sd,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)


g_artifact_sd <-
  graph_from_adjacency_matrix(sim_artifact_sd,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

g_rand_t_sd <-
  graph_from_adjacency_matrix(sim_rand_t_sd,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

g_rand_b_sd <-
  graph_from_adjacency_matrix(sim_rand_b_sd,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)


# Project with Jaccard ----------------------------------------------------

g_prov_jacc <-
  graph_from_adjacency_matrix(sim_prov_jacc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)


g_artifact_jacc <-
  graph_from_adjacency_matrix(
    sim_artifact_jacc,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )


g_rand_t_jacc <-
  graph_from_adjacency_matrix(sim_rand_t_jacc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

g_rand_b_jacc <-
  graph_from_adjacency_matrix(sim_rand_b_jacc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)
