# ============================================================
# Simulation: Iterative Density Filtering + DualSimplex
# ============================================================
# Description:
#   Reproducible simulation pipeline demonstrating the effect
#   of iterative density-based filtering on reference-free
#   deconvolution under noisy conditions.
#
# Steps:
#   1. Generate simulated bulk expression data
#   2. Add noise
#   3. Perform projection and density annotation
#   4. Apply iterative density filtering
#   5. Run DualSimplex optimization
#   6. Evaluate deconvolution accuracy
# ============================================================

set.seed(3)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
n_ct <- 3
n_genes <- 10000
n_samples <- 100

# ------------------------------------------------------------
# Simulation
# ------------------------------------------------------------
sim <- create_simulation(
  n_genes = n_genes,
  n_samples = n_samples,
  n_cell_types = n_ct,
  with_marker_genes = FALSE
)

true_basis <- sim$basis
true_proportions <- sim$proportions

# ------------------------------------------------------------
# Noise injection
# ------------------------------------------------------------
sim <- sim %>% add_noise(noise_deviation = 5)
data_raw <- sim$data

# ------------------------------------------------------------
# DualSimplex initialization
# ------------------------------------------------------------
dso <- DualSimplexSolver$new()
dso$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
dso$project(n_ct)

dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

# ------------------------------------------------------------
# Density annotation (before filtering)
# ------------------------------------------------------------
add_annov2(n_ct = n_ct, genes = TRUE)

density_values <- dso[["st"]][["data"]]@featureData@data[["density"]]
df <- data.frame(density = density_values)

ggplot(df, aes(x = density)) +
  geom_histogram(bins = 200, fill = "lightblue", color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Gene density distribution (before filtering)",
    x = "Density",
    y = "Frequency"
  )

# ------------------------------------------------------------
# Iterative density filtering
# ------------------------------------------------------------
dso <- iterative_density_filter(
  dso,
  threshold = 1,
  n_cell_types = n_ct,
  genes = TRUE
)

# ------------------------------------------------------------
# Projection after filtering
# ------------------------------------------------------------
dso$project(n_ct)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

# ------------------------------------------------------------
# Optimization and evaluation
# ------------------------------------------------------------
dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()

ptp <- coerce_pred_true_props(
  solution$H,
  true_proportions
)

plot_correlation_matrix(ptp[[1]], ptp[[2]])
plot_ptp_lines(ptp)
plot_ptp_scatter(ptp)
