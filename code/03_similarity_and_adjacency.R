
# Run previous scripts ----------------------------------------------------

source("code/01_bipartite_graph.R")
source("code/02_similarity_functions.R")

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

