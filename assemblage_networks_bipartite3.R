
require(tidyverse)
require(igraph)
require(ggraph)

# Import data from file ---------------------------------------------------

dat <- read_csv("Catalog_SiteA.csv",
                col_select = c(LEVEL_ID, CODE))

# Create un-weighted bipartite graph --------------------------------------

g_assemblages_bpg <-
  graph_from_data_frame(unique.data.frame(select(dat, LEVEL_ID, CODE)),
                        directed = FALSE)

V(g_assemblages_bpg)$type <-
  bipartite_mapping(g_assemblages_bpg)$type

g_assemblages_bpg_layout <-
  g_assemblages_bpg %>% layout_as_bipartite()

g_assemblages_bpg %>%
  ggraph(layout = g_assemblages_bpg_layout) +
  geom_edge_link(color = "gray", alpha = 0.25) +
  geom_node_point(aes(color = type)) +
  scale_color_manual(
    values = c("green", "blue"),
    name = "Node Type",
    breaks = c(FALSE, TRUE),
    labels = c("Provenience", "Artifact")
  ) +
  ggtitle("Bipartite network of Provenience and Artifact Type")


# Create incidence matrix from bipartite graph -----------------------------

g_assemblages_bpg_inc <- as_incidence_matrix(g_assemblages_bpg)

# Project one-mode graphs with igraph --------------------------------------

g_assemblages_proj <-
  bipartite_projection(g_assemblages_bpg, multiplicity = TRUE)

g_assemblages_proj_prov <- g_assemblages_proj$proj1

g_assemblages_proj_artifact <- g_assemblages_proj$proj2

g_assemblages_proj_prov %>%
  ggraph() +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "green") +
  ggtitle("Network of Proveniences")

g_assemblages_proj_artifact %>%
  ggraph() +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "blue") +
  ggtitle("Network of Artifact Types")


# Project one-mode graphs, Szymkiewicz-Simpson -----------------------------

overlap_coef_bin <- function(x) {
  # Calculate the pairwise sums of non-zero matrix elements to find the number
  # of intersecting elements between each input column
  bin_intersect_mat <- t(x) %*% x
  
  # Calculate the input column sums to find the individual size of each set
  col_sum <- apply(x, 2, function(xx)
    sum(xx != 0))
  
  # Find the smaller of each pair of sets by taking the matrix outer minimum
  min_set_size_mat <- outer(col_sum, col_sum, FUN = pmin)
  
  # Szymkiewicz-Simpson overlap coefficient is the pairwise intersection of two
  # sets divided by the size of the smaller set
  res <- bin_intersect_mat / min_set_size_mat
  
  # Set diagonal to identity
  diag(res) <- 1L
  
  # Assign input column names to rows and columns of return matrix
  dimnames(res) <- list(colnames(x), colnames(x))
  
  return(res)
}


## Project provenience -----------------------------------------------------

prov_sim_ssoc <- overlap_coef_bin(t(g_assemblages_bpg_inc))

prov_ssoc_sims <-
  prov_sim_ssoc[lower.tri(prov_sim_ssoc, diag = FALSE)]

ggplot(data = data.frame(x = c(prov_ssoc_sims)), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Szymkiewicz-Simpson Similarity for Proveniences")

g_assemblages_proj_prov_oc <-
  graph_from_adjacency_matrix(prov_sim_ssoc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

g_assemblages_proj_prov_oc %>%
  ggraph(layout = "mds") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "green") +
  ggtitle("Network of Proveniences")

ggplot(data = data.frame(x = degree(g_assemblages_proj_prov_oc)), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Szymkiewicz-Simpson Degree for Proveniences")

ggplot(data = data.frame(x = strength(g_assemblages_proj_prov_oc)), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Szymkiewicz-Simpson Edge-weighted Degree for Proveniences")

ggplot(data = data.frame(x = E(g_assemblages_proj_prov_oc)$weight), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Szymkiewicz-Simpson Edge Weight for Proveniences")


## Project artifact types --------------------------------------------------

artifact_sim_ssoc <- overlap_coef_bin(g_assemblages_bpg_inc)

artifact_ssoc_sims <-
  artifact_sim_ssoc[lower.tri(artifact_sim_ssoc, diag = FALSE)]

ggplot(data = data.frame(x = c(artifact_ssoc_sims)), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Szymkiewicz-Simpson Similarity for Artifacts")

g_assemblages_proj_artifact_oc <-
  graph_from_adjacency_matrix(
    artifact_sim_ssoc^3,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

g_assemblages_proj_artifact_oc %>%
  ggraph(layout = "kk") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "blue") +
  ggtitle("Network of Artifacts")

ggplot(data = data.frame(x = degree(g_assemblages_proj_artifact_oc)), aes(x = x)) +
  geom_density(fill = "blue") +
  ggtitle("Distribution of Szymkiewicz-Simpson Degree for Artifacts")

ggplot(data = data.frame(x = strength(g_assemblages_proj_artifact_oc)), aes(x = x)) +
  geom_density(fill = "blue") +
  ggtitle("Distribution of Szymkiewicz-Simpson Node Strength for Artifacts")

ggplot(data = data.frame(x = E(g_assemblages_proj_artifact_oc)$weight), aes(x = x)) +
  geom_density(fill = "blue") +
  ggtitle("Distribution of Szymkiewicz-Simpson Edge Weight for Artifacts")


# Project one-mode graphs, Sorenson-Dice -----------------------------------

soren_dice_sim_bin <- function(x) {
  # Calculate the pairwise sums of non-zero matrix elements to find the number
  # of intersecting elements between each input column
  bin_intersect_mat <- t(x) %*% x
  
  # Calculate the input column sums, to find the individual size of each set
  col_sum <- apply(x, 2, function(xx)
    sum(xx != 0))
  
  # Calculate the matrix outer sums for pairwise sum of set sizes
  set_size_sum_mat <- outer(col_sum, col_sum, FUN = "+")
  
  # Sorenson-Dice index is twice the size of the intersection divided by the
  # sum of the size for each set
  res <- (2 * bin_intersect_mat) / set_size_sum_mat
  
  # Set diagonal to identity
  diag(res) <- 1L
  
  # Assign input column names to rows and columns of return matrix
  dimnames(res) <- list(colnames(x), colnames(x))
  return(res)
  
}

## Project provenience -----------------------------------------------------

prov_sim_sd <- soren_dice_sim_bin(t(g_assemblages_bpg_inc))

prov_sd_sims <-
  prov_sim_sd[lower.tri(prov_sim_sd, diag = FALSE)]

ggplot(data = data.frame(x = prov_sd_sims), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Sorenson-Dice Similarity for Provenience")

g_assemblages_proj_prov_sd <-
  graph_from_adjacency_matrix(prov_sim_sd,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

# g_assemblages_proj_prov_sd <-
#   graph_from_adjacency_matrix(ifelse(prov_adj_sd < 0.439, 0, prov_adj_sd),
#                               mode = "undirected",
#                               weighted = TRUE,
#                               diag = FALSE)

# g_assemblages_proj_prov_sd %>%
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "green") +
#   ggtitle("Network of Proveniences")

data.frame(x = degree(g_assemblages_proj_prov_sd)) %>% filter(x > 0) %>%
  ggplot(aes(x = x)) +
  geom_histogram(color = "gray",
                 fill = "green",
                 bins = 20)

ggplot(data = data.frame(x = degree(g_assemblages_proj_prov_sd)), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Sorenson-Dice Degree for Proveniences")

ggplot(data = data.frame(x = E(g_assemblages_proj_prov_sd)$weight), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Sorenson-Dice Edge Weight for Proveniences")


## Project artifact types --------------------------------------------------

artifact_sim_sd <- soren_dice_sim_bin(g_assemblages_bpg_inc)

artifact_sd_sims <-
  artifact_sim_sd[lower.tri(artifact_sim_sd, diag = FALSE)]

ggplot(data = data.frame(x = artifact_sd_sims), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Sorenson-Dice Similarity for Artifacts")

g_assemblages_proj_artifact_sd <-
  graph_from_adjacency_matrix(artifact_sim_sd,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

# g_assemblages_artifact_prov_sd %>%
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "green") +
#   ggtitle("Network of Artifacts")

ggplot(data = data.frame(x = degree(g_assemblages_proj_artifact_sd)), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Sorenson-Dice Degree for Artifacts")

data.frame(x = degree(g_assemblages_proj_artifact_sd)) %>%
  ggplot(aes(x = x)) +
  geom_histogram(color = "gray",
                 fill = "blue",
                 bins = 20)


ggplot(data = data.frame(x = E(g_assemblages_proj_artifact_sd)$weight), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Sorenson-Dice Edge Weight for Artifacts")


# Project one-mode graphs, Jaccard ----------------------------------------

jaccard_sim_bin <- function(x) {
  # Calculate the pairwise sums of non-zero matrix elements to find the number
  # of intersecting elements between each input column
  bin_intersect_mat <- t(x) %*% x
  
  # Calculate the input column sums to find the individual size of each set
  x_col_sum <- apply(x, 2, function(xx)
    sum(xx != 0))
  
  # Calculate the matrix outer sums for pairwise sum of set sizes
  set_size_sum_mat <- outer(x_col_sum, x_col_sum, FUN = "+")
  
  # Jaccard index is intersection of set sizes over the size of the union of
  # sets
  res <- bin_intersect_mat / (set_size_sum_mat - bin_intersect_mat)
  
  # Set diagonal to identity
  diag(res) <- 1L
  
  # Assign input column names to rows and columns of return matrix
  dimnames(res) <- list(colnames(x), colnames(x))
  
  return(res)
}

## Project provenience -----------------------------------------------------

prov_sim_jacc <- jaccard_sim_bin(t(g_assemblages_bpg_inc))

prov_jacc_sims <-
  prov_sim_jacc[lower.tri(prov_sim_jacc, diag = FALSE)]

ggplot(data = data.frame(x = prov_jacc_sims), aes(x = x)) +
  geom_density(fill = "green") +
  ggtitle("Distribution of Jaccard Similarity for Provenience")

g_assemblages_proj_prov_jacc <-
  graph_from_adjacency_matrix(prov_sim_jacc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

# g_assemblages_proj_prov_jacc %>%
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "green") +
#   ggtitle("Network of Proveniences")

ggplot(data = data.frame(x = degree(g_assemblages_proj_prov_jacc)), aes(x = x)) +
  geom_density(fill = "green") +
  ggtitle("Distribution of Jaccard Degree for Proveniences")

ggplot(data = data.frame(x = E(g_assemblages_proj_prov_jacc)$weight), aes(x = x)) +
  geom_density(fill = "green") +
  ggtitle("Distribution of Jaccard Edge Weight for Proveniences")


## Project artifact types --------------------------------------------------

artifact_sim_jacc <- jaccard_sim_bin(g_assemblages_bpg_inc)

artifact_jacc_sims <-
  artifact_sim_jacc[lower.tri(artifact_sim_jacc, diag = FALSE)]

ggplot(data = data.frame(x = artifact_jacc_sims), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Jaccard Similarity for Artifacts")

g_assemblages_proj_artifact_jacc <-
  graph_from_adjacency_matrix(
    artifact_sim_jacc,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

# g_assemblages_proj_prov_jacc %>%
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "blue") +
#   ggtitle("Network of Artifacts")

ggplot(data = data.frame(x = degree(g_assemblages_proj_artifact_jacc)), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Jaccard Degree for Artifacts")

ggplot(data = data.frame(x = E(g_assemblages_proj_artifact_jacc)$weight), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Jaccard Edge Weight for Artifacts")

# TOM Adjacency -----------------------------------------------------------

tom_similarity <- function(adj_mat) {
  l_mat <- adj_mat %*% t(adj_mat)
  
  k_row <- rowSums(adj_mat)
  k_col <- colSums(adj_mat)
  k_min <- outer(k_row, k_col, FUN = pmin)
  
  tom_sim <- (l_mat + adj_mat) / (k_min + 1 - adj_mat)
  
  diag(tom_sim) <- 1
  
  dimnames(tom_sim) <- list(colnames(adj_mat), colnames(adj_mat))
  
  return(tom_sim)
}

## Provenience overlap matrix ----------------------------------------------



## Artifact type overlap matrix --------------------------------------------


# Distribution of sample sizes ---------------------------------------------

require(gridExtra)

grid.arrange(
  ggplot(data = data.frame(x = rowSums(
    g_assemblages_bpg_inc
  )),
  aes(x = x, y = after_stat(density))) +
    geom_histogram(
      color = "gray",
      fill = "green",
      bins = 20
    ) +
    ggtitle("Artifact Types per Provenience"),
  
  ggplot(data = data.frame(x = colSums(
    g_assemblages_bpg_inc
  )),
  aes(x = x, y = after_stat(density))) +
    geom_histogram(
      color = "gray",
      fill = "blue",
      bins = 20
    ) +
    ggtitle("Occurence per Artifact Type"),
  
  ncol = 2
)



# Comparing similarity measures -------------------------------------------

prov_sims <-
  data.frame(ssoc = prov_ssoc_sims,
             jacc = prov_jacc_sims,
             sd = prov_sd_sims)

prov_sims %>% stack() %>%
  # filter(values > 0) %>% # View non-zero entries
  ggplot(aes(x = values)) +
  geom_density(color = "green",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "free")

prov_sims %>% stack() %>%
  filter(values > 0) %>%
  ggplot(aes(x = ind, y = values, fill = ind)) +
  geom_violin()

artifact_sims <-
  data.frame(ssoc = artifact_ssoc_sims,
             jacc = artifact_jacc_sims,
             sd = artifact_sd_sims)

artifact_sims %>% stack() %>%
  # filter(values > 0) %>% # View non-zero entries
  ggplot(aes(x = values)) +
  geom_density(color = "blue",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "free")

artifact_sims %>% stack() %>%
  # filter(values > 0) %>%
  ggplot(aes(x = ind, y = values, fill = ind)) +
  geom_violin()



# Thresholding ------------------------------------------------------------



## Hard threshold ---------------------------------------------------------

signum_adj <- function(x, tau = 0.5) {
  ifelse(x < tau, 0, 1)
}

## Soft threshold ---------------------------------------------------------

sigmoid_adj <- function(x, alpha = 10, tau = 0.5) {
  1.0 / (1.0 + exp(-alpha * (x - tau)))
}

power_adj <- function(x, beta = 1) {
  x ^ beta
}


# Scale-free topology -----------------------------------------------------



## Power law probability distribution -------------------------------------



powerlaw_pdf <- function(x, gamma = 1) {
  x ^ (-gamma)
}

ggplot(data.frame(x = c(1, 251)), aes(x = x)) +
  stat_function(fun = powerlaw_pdf, args = list(gamma = 1), aes(color = "1")) +
  stat_function(fun = powerlaw_pdf, args = list(gamma = 2), aes(color = "2")) +
  stat_function(fun = powerlaw_pdf, args = list(gamma = 3), aes(color = "3")) +
  stat_function(fun = powerlaw_pdf, args = list(gamma = 4), aes(color = "4")) +
  stat_function(fun = powerlaw_pdf, args = list(gamma = 5), aes(color = "5")) +
  scale_y_continuous(trans = "log10") +
  labs(color = "Gamma value") +
  ylab("log P(x)")


# Community detection -----------------------------------------------------

artifact_adj_power <- power_adj(artifact_sim_ssoc, beta = 3)

g_artifact <-
  graph_from_adjacency_matrix(
    artifact_adj_power,
    mode = 'undirected',
    weighted = TRUE,
    diag = FALSE
  )

prov_adj_signum <- signum_adj(prov_sim_jacc, tau = 0.40)

g_prov <-
  graph_from_adjacency_matrix(prov_adj_signum,
                              mode = 'undirected',
                              weighted = NULL,
                              diag = FALSE)

