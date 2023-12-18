source("code/00_data_import.R", echo = TRUE)

require(igraph)
require(WGCNA)

WGCNA::enableWGCNAThreads()


# Create un-weighted bipartite assemblage network -------------------------

g_assemblages_bpg <-
  graph_from_data_frame(unique.data.frame(select(dat, LEVEL_ID, CODE)),
                        directed = FALSE)


# Assign vertex types -----------------------------------------------------

V(g_assemblages_bpg)$type <-
  bipartite_mapping(g_assemblages_bpg)$type


# Create incidence matrix from bipartite graph -----------------------------

g_assemblages_bpg_inc <- as_biadjacency_matrix(g_assemblages_bpg)


# Create random bipartite graph (for comparisons) -------------------------

g_random_bpg <-
  sample_bipartite(
    n1 = 152,
    n2 = 233,
    type = "gnm",
    m = 2240,
    directed = FALSE
  )

g_random_bpg_inc <- as_biadjacency_matrix(g_random_bpg)

