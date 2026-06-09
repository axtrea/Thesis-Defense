# ============================================================
# Benchmark: DualSimplex with vs without filtering
# Noise robustness evaluation
# ============================================================

library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Noise wrappers (extendable design)
# ------------------------------------------------------------

add_poisson_noise_wrapper <- function(sim, noise_level) {
  sim <- add_poisson_noise(sim, lambda_noise = noise_level)
  return(sim)
}

add_gaussian_noise_wrapper <- function(sim, noise_level) {
  sim <- add_noise(sim, noise_deviation = noise_level)
  return(sim)
}

add_noise_factory <- function(noise_type) {
  switch(
    noise_type,
    poisson = add_poisson_noise_wrapper,
    gaussian = add_gaussian_noise_wrapper,
    stop("Unknown noise type")
  )
}

# ------------------------------------------------------------
# Core experiment function
# ------------------------------------------------------------

run_single_experiment <- function(
    noise_level,
    noise_type = "poisson",
    n_genes = 10000,
    n_samples = 100,
    n_cell_types = 3,
    seed = 1
) {

  set.seed(seed)

  # Create simulation
  sim <- create_simulation(
    n_genes = n_genes,
    n_samples = n_samples,
    n_cell_types = n_cell_types,
    with_marker_genes = FALSE
  )

  # Add noise (generic)
  noise_fn <- add_noise_factory(noise_type)
  sim <- noise_fn(sim, noise_level)

  data_raw <- sim$data
  true_props <- sim$proportions

  # ============================================================
  # WITHOUT FILTERING
  # ============================================================
  dso_raw <- DualSimplexSolver$new()
  dso_raw$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
  dso_raw$project(n_cell_types)
  dso_raw$init_solution("random")
  dso_raw$default_optimization()
  dso_raw$finalize_solution()

  solution_raw <- dso_raw$get_solution()
  ptp_raw <- coerce_pred_true_props(solution_raw$H, true_props)
  rmse_raw <- rmse_loss_function(ptp_raw[[1]], ptp_raw[[2]])

  # ============================================================
  # WITH FILTERING
  # ============================================================
  dso_filt <- DualSimplexSolver$new()
  dso_filt$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
  dso_filt$project(n_cell_types)

  dso_filt <- iterative_density_filter(
    dso_filt,
    threshold = 100,
    n_cell_types = n_cell_types,
    genes = TRUE,
    max_iter = 30,
    verbose = FALSE
  )

  dso_filt$init_solution("random")
  dso_filt$default_optimization()
  dso_filt$finalize_solution()

  solution_filt <- dso_filt$get_solution()
  ptp_filt <- coerce_pred_true_props(solution_filt$H, true_props)
  rmse_filt <- rmse_loss_function(ptp_filt[[1]], ptp_filt[[2]])

  # Return results
  data.frame(
    noise = noise_level,
    noise_type = noise_type,
    seed = seed,
    method = c("without filtering", "filtered"),
    rmse = c(rmse_raw, rmse_filt)
  )
}

# ------------------------------------------------------------
# Benchmark runner
# ------------------------------------------------------------

run_benchmark <- function(
    noise_levels,
    noise_type = "poisson",
    n_replicates = 5
) {

  results_list <- list()

  for (nl in noise_levels) {
    cat("\n=== Noise level:", nl, "| Type:", noise_type, "===\n")

    for (rep in 1:n_replicates) {

      res <- run_single_experiment(
        noise_level = nl,
        noise_type = noise_type,
        seed = rep
      )

      results_list[[length(results_list) + 1]] <- res
    }
  }

  bind_rows(results_list)
}

# ------------------------------------------------------------
# Run experiments
# ------------------------------------------------------------

noise_levels <- c(2, 4, 6, 8, 10)

results_all <- run_benchmark(
  noise_levels = noise_levels,
  noise_type = "poisson",
  n_replicates = 10
)

# ------------------------------------------------------------
# Visualization
# ------------------------------------------------------------

p <- ggplot(results_all, aes(x = as.factor(noise), y = rmse, fill = method)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_dodge(0.8), size = 1.1, alpha = 0.4, color = "gray30") +
  labs(
    title = paste("Noise robustness:", unique(results_all$noise_type)),
    x = "Noise level",
    y = "RMSE",
    fill = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Set2")

print(p)

# ------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------

summary_table <- results_all %>%
  group_by(noise, noise_type, method) %>%
  summarise(
    mean_rmse = mean(rmse),
    sd_rmse = sd(rmse),
    .groups = "drop"
  )

print(summary_table)
