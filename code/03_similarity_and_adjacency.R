


# Run previous scripts ----------------------------------------------------

source("code/01_bipartite_graph.R", echo = TRUE)
source("code/02_similarity_functions.R", echo = TRUE)

# Overlap coefficient similarity ------------------------------------------

sim_prov_oc <- overlap_coef_bin(t(g_assemblages_bpg_inc))

prov_ssoc_sims <-
  sim_prov_oc[lower.tri(sim_prov_oc, diag = FALSE)]

sim_artifact_oc <- overlap_coef_bin(g_assemblages_bpg_inc)

artifact_ssoc_sims <-
  sim_artifact_oc[lower.tri(sim_artifact_oc, diag = FALSE)]

sim_rand_t_oc <- overlap_coef_bin(t(g_random_bpg_inc))

rand_t_ssoc_sims <-
  sim_rand_t_oc[lower.tri(sim_rand_t_oc, diag = FALSE)]

sim_rand_b_oc <- overlap_coef_bin(g_random_bpg_inc)

rand_b_ssoc_sims <-
  sim_rand_b_oc[lower.tri(sim_rand_b_oc, diag = FALSE)]



# Sorenson-Dice similarity ------------------------------------------------


sim_prov_sd <- soren_dice_sim_bin(t(g_assemblages_bpg_inc))

prov_sd_sims <-
  sim_prov_sd[lower.tri(sim_prov_sd, diag = FALSE)]

sim_artifact_sd <- soren_dice_sim_bin(g_assemblages_bpg_inc)

artifact_sd_sims <-
  sim_artifact_sd[lower.tri(sim_artifact_sd, diag = FALSE)]

sim_rand_t_sd <- soren_dice_sim_bin(t(g_random_bpg_inc))

rand_t_sd_sims <-
  sim_rand_t_sd[lower.tri(sim_rand_t_sd, diag = FALSE)]

sim_rand_b_sd <- soren_dice_sim_bin(g_random_bpg_inc)

rand_b_sd_sims <-
  sim_rand_b_sd[lower.tri(sim_rand_b_sd, diag = FALSE)]


# Jaccard similarity ------------------------------------------------------



sim_prov_jacc <- jaccard_sim_bin(t(g_assemblages_bpg_inc))

prov_jacc_sims <-
  sim_prov_jacc[lower.tri(sim_prov_jacc, diag = FALSE)]


sim_artifact_jacc <- jaccard_sim_bin(g_assemblages_bpg_inc)

artifact_jacc_sims <-
  sim_artifact_jacc[lower.tri(sim_artifact_jacc, diag = FALSE)]


sim_rand_t_jacc <- jaccard_sim_bin(t(g_random_bpg_inc))

rand_t_jacc_sims <-
  sim_rand_t_jacc[lower.tri(sim_rand_t_jacc, diag = FALSE)]


sim_rand_b_jacc <- jaccard_sim_bin(g_random_bpg_inc)

rand_b_jacc_sims <-
  sim_rand_b_jacc[lower.tri(sim_rand_b_jacc, diag = FALSE)]


# Comparing similarity measures -------------------------------------------


## Archaeological graphs --------------------------------------------------


prov_sims <-
  data.frame(ssoc = prov_ssoc_sims,
             jacc = prov_jacc_sims,
             sd = prov_sd_sims)

prov_sims %>% stack() %>%
  # filter(values > 0) %>% # View non-zero entries
  ggplot(aes(x = values)) +
  geom_density(fill = "darkgreen",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Similarity for Provenience")

prov_sims %>% stack() %>%
  # filter(values > 0) %>%
  ggplot(aes(x = ind, y = values, fill = ind)) +
  geom_violin() +
  ggtitle("Distribution of Similarity for Provenience")

artifact_sims <-
  data.frame(ssoc = artifact_ssoc_sims,
             jacc = artifact_jacc_sims,
             sd = artifact_sd_sims)

artifact_sims %>% stack() %>%
  # filter(values > 0) %>% # View non-zero entries
  ggplot(aes(x = values)) +
  geom_density(fill = "darkblue",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Similarity for Artifacts")

artifact_sims %>% stack() %>%
  # filter(values > 0) %>%
  ggplot(aes(x = ind, y = values, fill = ind)) +
  geom_violin() +
  ggtitle("Distribution of Similarity for Artifacts")


## Random graph distributions----------------------------------------------

rand_t_sims <-
  data.frame(ssoc = rand_t_ssoc_sims,
             jacc = rand_t_jacc_sims,
             sd = rand_t_sd_sims)

rand_t_sims %>% stack() %>%
  # filter(values > 0) %>% # View non-zero entries
  ggplot(aes(x = values)) +
  geom_density(fill = "darkgreen",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Similarity for Random 'Provenience'")

rand_t_sims %>% stack() %>%
  # filter(values > 0) %>%
  ggplot(aes(x = ind, y = values, fill = ind)) +
  geom_violin() +
  ggtitle("Distribution of Similarity for Random 'Provenience'")

rand_b_sims <-
  data.frame(ssoc = rand_b_ssoc_sims,
             jacc = rand_b_jacc_sims,
             sd = rand_b_sd_sims)

rand_b_sims %>% stack() %>%
  # filter(values > 0) %>% # View non-zero entries
  ggplot(aes(x = values)) +
  geom_density(fill = "darkblue",
               alpha = 0.4) +
  facet_grid(ind ~ ., scales = "fixed") +
  ggtitle("Distribution of Similarity for Random 'Artifacts'")

rand_b_sims %>% stack() %>%
  # filter(values > 0) %>%
  ggplot(aes(x = ind, y = values, fill = ind)) +
  geom_violin() +
  ggtitle("Distribution of Similarity for Random 'Artifacts'")
