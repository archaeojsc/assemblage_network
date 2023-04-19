
# Required libraries ------------------------------------------------------


require(tidyverse)
require(igraph)
require(ggraph)

# Import data from file ---------------------------------------------------

# Data for Site A
dat_raw <- read_csv("./data/Catalog_SiteA.csv",
                col_select = c(LEVEL_ID, CODE))

# Data for Sites A and C
# dat_raw <- read_csv("./data/Catalog_AC.csv",
#                     col_select = c(LEVEL_ID, CODE))

# Artifacts to exclude from analysis
exclude_artifacts <-
  as.vector(read.csv("./data/code_exclude.csv", header = T)$x)

# Remove Bakelite and plastic buttons from exclusion list
exclude_artifacts <-
  exclude_artifacts[!exclude_artifacts %in% c("BKLT", "PB")]

# Import list of artifact codes
artifact_codes <-
  read.csv("./data/code_list.csv",
           header = TRUE,
           stringsAsFactors = FALSE)

# Filter data set for excluded artifact types
dat <- dat_raw %>% filter(!(CODE %in% exclude_artifacts))


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
  geom_edge_link(color = "darkgray", alpha = 0.25) +
  geom_node_point(aes(color = type)) +
  scale_color_manual(
    values = c("darkgreen", "darkblue"),
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
  ggraph(layout = "kk") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "darkgreen") +
  ggtitle("Network of Proveniences")

g_assemblages_proj_artifact %>%
  ggraph(layout = "kk") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "darkblue") +
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

sim_oc_prov <- overlap_coef_bin(t(g_assemblages_bpg_inc))

prov_ssoc_sims <-
  sim_oc_prov[lower.tri(sim_oc_prov, diag = FALSE)]

ggplot(data = data.frame(x = c(prov_ssoc_sims)), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Szymkiewicz-Simpson Similarity for Proveniences")

g_assemblages_proj_prov_oc <-
  graph_from_adjacency_matrix(sim_oc_prov,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

g_assemblages_proj_prov_oc %>%
  ggraph(layout = "auto") +
  geom_edge_link(color = "darkgray", aes(alpha = weight)) +
  geom_node_point(color = "darkgreen") +
  ggtitle("Network of Proveniences")

ggplot(data = data.frame(x = degree(g_assemblages_proj_prov_oc)), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Szymkiewicz-Simpson Degree for Proveniences")

ggplot(data = data.frame(x = strength(g_assemblages_proj_prov_oc)), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Szymkiewicz-Simpson Node Strength for Proveniences")

ggplot(data = data.frame(x = E(g_assemblages_proj_prov_oc)$weight), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Szymkiewicz-Simpson Edge Weight for Proveniences")


## Project artifact types --------------------------------------------------

sim_oc_artifact <- overlap_coef_bin(g_assemblages_bpg_inc)

artifact_ssoc_sims <-
  sim_oc_artifact[lower.tri(sim_oc_artifact, diag = FALSE)]

ggplot(data = data.frame(x = c(artifact_ssoc_sims^3)), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
  ggtitle("Distribution of Szymkiewicz-Simpson Similarity for Artifacts")

g_assemblages_proj_artifact_oc <-
  graph_from_adjacency_matrix(
    sim_oc_artifact,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

g_assemblages_proj_artifact_oc %>%
  ggraph(layout = "auto") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "darkblue") +
  ggtitle("Network of Artifacts")

ggplot(data = data.frame(x = degree(g_assemblages_proj_artifact_oc)), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
  ggtitle("Distribution of Szymkiewicz-Simpson Degree for Artifacts")

ggplot(data = data.frame(x = strength(g_assemblages_proj_artifact_oc)), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
  ggtitle("Distribution of Szymkiewicz-Simpson Node Strength for Artifacts")

ggplot(data = data.frame(x = E(g_assemblages_proj_artifact_oc)$weight), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
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

sim_sd_prov <- soren_dice_sim_bin(t(g_assemblages_bpg_inc))

prov_sd_sims <-
  sim_sd_prov[lower.tri(sim_sd_prov, diag = FALSE)]

ggplot(data = data.frame(x = prov_sd_sims), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Sorenson-Dice Similarity for Provenience")

g_assemblages_proj_prov_sd <-
  graph_from_adjacency_matrix(sim_sd_prov,
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
#   geom_edge_link(color = "darkgray", aes(alpha = weight)) +
#   geom_node_point(color = "darkgreen") +
#   ggtitle("Network of Proveniences")

data.frame(x = degree(g_assemblages_proj_prov_sd)) %>% 
  ggplot(aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4)+
  ggtitle("Distribution of Sorenson-Dice Node Degree for Proveniences")

ggplot(data = data.frame(x = strength(g_assemblages_proj_prov_sd)), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Sorenson-Dice Node Strength for Proveniences")

ggplot(data = data.frame(x = E(g_assemblages_proj_prov_sd)$weight), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Sorenson-Dice Edge Weight for Proveniences")


## Project artifact types --------------------------------------------------

sim_sd_artifact <- soren_dice_sim_bin(g_assemblages_bpg_inc)

artifact_sd_sims <-
  sim_sd_artifact[lower.tri(sim_sd_artifact, diag = FALSE)]

ggplot(data = data.frame(x = artifact_sd_sims), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
  ggtitle("Distribution of Sorenson-Dice Similarity for Artifacts")

g_assemblages_proj_artifact_sd <-
  graph_from_adjacency_matrix(sim_sd_artifact,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

# g_assemblages_artifact_prov_sd %>%
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "darkgray", aes(alpha = weight)) +
#   geom_node_point(color = "darkgreen") +
#   ggtitle("Network of Artifacts")

ggplot(data = data.frame(x = degree(g_assemblages_proj_artifact_sd)), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
  ggtitle("Distribution of Sorenson-Dice Degree for Artifacts")

data.frame(x = strength(g_assemblages_proj_artifact_sd)) %>%
  ggplot(aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4)+
  ggtitle("Distribution of Sorenson-Dice Node Strength for Artifacts")


ggplot(data = data.frame(x = E(g_assemblages_proj_artifact_sd)$weight), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
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

sim_jacc_prov <- jaccard_sim_bin(t(g_assemblages_bpg_inc))

prov_jacc_sims <-
  sim_jacc_prov[lower.tri(sim_jacc_prov, diag = FALSE)]

ggplot(data = data.frame(x = prov_jacc_sims), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Jaccard Similarity for Provenience")

g_assemblages_proj_prov_jacc <-
  graph_from_adjacency_matrix(sim_jacc_prov,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

# g_assemblages_proj_prov_jacc %>%
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "darkgray", aes(alpha = weight)) +
#   geom_node_point(color = "darkgreen") +
#   ggtitle("Network of Proveniences")

ggplot(data = data.frame(x = degree(g_assemblages_proj_prov_jacc)), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Jaccard Degree for Proveniences")

ggplot(data = data.frame(x = strength(g_assemblages_proj_prov_jacc)), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Jaccard Node Strength for Proveniences")

ggplot(data = data.frame(x = E(g_assemblages_proj_prov_jacc)$weight), aes(x = x)) +
  geom_density(fill = "darkgreen", alpha = 0.4) +
  ggtitle("Distribution of Jaccard Edge Weight for Proveniences")


## Project artifact types --------------------------------------------------

sim_jacc_artifact <- jaccard_sim_bin(g_assemblages_bpg_inc)

artifact_jacc_sims <-
  sim_jacc_artifact[lower.tri(sim_jacc_artifact, diag = FALSE)]

ggplot(data = data.frame(x = artifact_jacc_sims), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
  ggtitle("Distribution of Jaccard Similarity for Artifacts")

g_assemblages_proj_artifact_jacc <-
  graph_from_adjacency_matrix(
    sim_jacc_artifact,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

# g_assemblages_proj_prov_jacc %>%
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "darkgray", aes(alpha = weight)) +
#   geom_node_point(color = "darkblue") +
#   ggtitle("Network of Artifacts")

ggplot(data = data.frame(x = degree(g_assemblages_proj_artifact_jacc)), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
  ggtitle("Distribution of Jaccard Degree for Artifacts")

ggplot(data = data.frame(x = strength(g_assemblages_proj_artifact_jacc)), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
  ggtitle("Distribution of Jaccard Node Strength for Artifacts")

ggplot(data = data.frame(x = E(g_assemblages_proj_artifact_jacc)$weight), aes(x = x)) +
  geom_density(fill = "darkblue", alpha = 0.4) +
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
      color = "darkgray",
      fill = "darkgreen",
      bins = 20
    ) +
    ggtitle("Artifact Types per Provenience"),
  
  ggplot(data = data.frame(x = colSums(
    g_assemblages_bpg_inc
  )),
  aes(x = x, y = after_stat(density))) +
    geom_histogram(
      color = "darkgray",
      fill = "darkblue",
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
  geom_density(fill = "darkgreen",
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
  geom_density(fill = "darkblue",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "free")

artifact_sims %>% stack() %>%
  # filter(values > 0) %>%
  ggplot(aes(x = ind, y = values, fill = ind)) +
  geom_violin()



# Thresholding ------------------------------------------------------------
require(WGCNA)


## Hard threshold ---------------------------------------------------------

signum_adj <- function(x, tau = 0.5) {
  ifelse(x < tau, 0, 1)
}

## Soft threshold ---------------------------------------------------------


### Sigmoid adjacency -----------------------------------------------------

sigmoid_adj <- function(x, alpha, mu) {
  1.0 / (1.0 + exp(-alpha * (x - mu)))
}

#### Sigmoid adjacency example plot ---------------------------------------

ggplot(data.frame(x = seq(-10, 10, length.out = 1000)), aes(x = x)) +
  stat_function(fun = sigmoid_adj,
                args = list(alpha = 1, mu = 0),
                aes(color = "1|0")) +
  stat_function(fun = sigmoid_adj,
                args = list(alpha = 2, mu = 0),
                aes(color = "2|0")) +
  stat_function(fun = sigmoid_adj,
                args = list(alpha = 1, mu = 3),
                aes(color = "1|3")) +
  stat_function(fun = sigmoid_adj,
                args = list(alpha = 2, mu = 3),
                aes(color = "2|3")) +
  labs(color = "alpha|mu") +
  ylab("Sigmoid(x)") + xlab("x")



### Power adjacency --------------------------------------------------------

power_adj <- function(x, beta = 1) {
  abs(x) ^ beta
}


#### Power adjacency example plot -------------------------------------------

ggplot(data.frame(x = seq(0, 1, length.out = 1000)), aes(x = x)) +
  stat_function(fun = power_adj, args = list(beta = 2), aes(color = "2")) +
  stat_function(fun = power_adj, args = list(beta = 3), aes(color = "3")) +
  stat_function(fun = power_adj, args = list(beta = 4), aes(color = "4")) +
  stat_function(fun = power_adj, args = list(beta = 5), aes(color = "5")) +
  stat_function(fun = power_adj, args = list(beta = 6), aes(color = "6")) +
  labs(color = "beta") +
  ylab("Power(x)") + xlab("x")



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

artifact_adj <- power_adj(artifact_sim_ssoc, beta = 5)

# artifact_adj <- signum_adj(artifact_sim_ssoc, tau = 0.7)

# artifact_adj<- signum_adj(artifact_sim_sd, tau = 0.55)

g_artifact <-
  graph_from_adjacency_matrix(
    artifact_adj,
    mode = 'undirected',
    weighted = TRUE,
    diag = FALSE
  )

# Weighted graph 
g_artifact %>%
  ggraph(layout = "mds") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "darkblue")

# Unweighted graph 
# g_artifact %>%
#   ggraph(layout = "mds") +
#   geom_edge_link(color = "gray", alpha = 0.4) +
#   geom_node_point(color = "darkblue")


g_artifact2 <-
  graph_from_adjacency_matrix(
    signum_adj(artifact_sim_sd, tau = 0.55),
    mode = 'undirected',
    weighted = NULL,
    diag = FALSE
  )

g_artifact2 %>%
  ggraph(layout = "kk") +
  geom_edge_link(color = "gray", alpha = 0.4) +
  geom_node_point(color = "darkblue")


# prov_adj_signum <- power_adj(prov_sim_ssoc, beta = 2)

g_prov <-
  graph_from_adjacency_matrix(prov_adj_signum,
                              mode = 'undirected',
                              weighted = TRUE,
                              diag = FALSE)

g_prov %>%
  ggraph(layout = "auto") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "darkgreen")



# Testing scale-free adjacency construction -------------------------------

sim_mat_artifact <-
  c('artifact_sim_ssoc', 'artifact_sim_sd', 'artifact_sim_jacc')

sim_mat_prov <- c('prov_sim_ssoc', 'prov_sim_sd', 'prov_sim_jacc')

g_test <- g_assemblages_proj_artifact_oc

thresh.vals <- seq(0.6, 0.8, by = 0.001)
test.p<-c()

for (i in 1:length(thresh.vals)) {
  test.p[i] <-
    fit_power_law(degree(delete.edges(
      g_test, which(E(g_test)$weight < thresh.vals[i])
    )))$KS.stat
}

plot(thresh.vals, test.p)
abline(v=thresh.vals[which.min(test.p[test.p>0])], col = "red")


ggplot(data = data.frame(x = degree(delete.edges(
  g_test, which(E(g_test)$weight < thresh.vals[which.min(test.p[test.p>0])])
))), aes(x = x)) +
  geom_density()


# # testing Soft Threshold -------------------------------------------------
# 
# old.par <- par(no.readonly = TRUE)
# 
# test <-
#   pickSoftThreshold.fromSimilarity(prov_sim_ssoc,
#                                    verbose = 5,
#                                    moreNetworkConcepts = TRUE)
# 
# par(mfrow = c(1, 2))
# 
# plot(
#   test$fitIndices[, 1],
#   -sign(test$fitIndices[, 3]) * test$fitIndices[, 4],
#   xlab = "Soft Threshold (power)",
#   ylab = "Scale Free Topology Model Fit, signed R^2",
#   main = paste("Scale independence"),
#   type = 'n'
# )
# 
# text(
#   test$fitIndices[, 1],
#   -sign(test$fitIndices[, 3]) * test$fitIndices[, 4],
#   labels = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
#   cex = 0.9,
#   col = "red"
# )
# 
# abline(h = 0.85, col = "red")
# 
# plot(
#   test$fitIndices[, 1],
#   test$fitIndices[, 5],
#   xlab = "Soft Threshold (power)",
#   ylab = "Mean Connectivity",
#   type = "n",
#   main = paste("Mean connectivity")
# )
# 
# text(
#   test$fitIndices[, 1],
#   test$fitIndices[, 5],
#   labels = test$fitIndices[, 1],
#   cex = .9,
#   col = "red"
# )
# 
# par(old.par)
# 
# SoftPower <- 2
# 
# test.adj <-
#   adjacency.fromSimilarity(prov_sim_ssoc, power = SoftPower)
# 
# # Testing Hard Threshold --------------------------------------------------
# 
# old.par <- par(no.readonly = TRUE)
# 
# test <-
#   pickHardThreshold.fromSimilarity(prov_sim_ssoc,
#                                    moreNetworkConcepts = TRUE,
#                                    RsquaredCut = 0.80)
# 
# par(mfrow = c(1, 2))
# 
# plot(
#   test$fitIndices[, 1],
#   -sign(test$fitIndices[, 4]) * test$fitIndices[, 3],
#   xlab = "Hard Threshold",
#   ylab = "Scale Free Topology Model Fit, signed R^2",
#   main = paste("Scale independence"),
#   type = 'n'
# )
# 
# text(
#   test$fitIndices[, 1],
#   -sign(test$fitIndices[, 4]) * test$fitIndices[, 3],
#   labels = test$fitIndices[, 1],
#   cex = 0.9,
#   col = "red"
# )
# 
# abline(h = 0.85, col = "red")
# 
# plot(
#   test$fitIndices[, 1],
#   test$fitIndices[, 6],
#   xlab = "Hard Threshold",
#   ylab = "Mean Connectivity",
#   type = "n",
#   main = paste("Mean connectivity")
# )
# 
# text(
#   test$fitIndices[, 1],
#   test$fitIndices[, 6],
#   labels = test$fitIndices[, 1],
#   cex = .9,
#   col = "red"
# )
# 
# par(old.par)
# 
# sigNum <- 0.25
# 
# test.adj <-
#   signumAdjacencyFunction(prov_sim_ssoc, threshold = sigNum)
# 
# # Testing TOM clustering --------------------------------------------------
# 
# 
# 
# test.TOM <- TOMsimilarity(test.adj)
# test.dissTOM <- 1 - test.TOM
# 
# testTree <- hclust(as.dist(test.dissTOM), method = "average")
# 
# plot(
#   testTree,
#   xlab = "",
#   sub = "",
#   main = "Clustering on TOM-based dissimilarity",
#   labels = FALSE,
#   hang = 0.04
# )
# 
# # Module identification using dynamic tree cut:
# dynamicMods = cutreeDynamic(
#   dendro = testTree,
#   distM = test.dissTOM,
#   deepSplit = 2,
#   pamRespectsDendro = FALSE
# )
# 
# table(dynamicMods)
# 
# # Convert numeric labels into colors
# dynamicColors = labels2colors(dynamicMods)
# table(dynamicColors)
# 
# plotDendroAndColors(
#   testTree,
#   dynamicColors,
#   "Dynamic Tree Cut",
#   dendroLabels = FALSE,
#   hang = 0.03,
#   addGuide = TRUE,
#   guideHang = 0.05,
#   main = "Dendrogram and module colors"
# )
# 

