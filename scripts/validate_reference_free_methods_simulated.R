# ======================================================================
# Validation of reference-free deconvolution methods on simulated data
# ======================================================================
# Simulation:
#   Genes        : 10,000
#   Samples      : 100
#   Cell types   : 15
#   Noise model  : Gaussian (sd = 5)
#
# Methods:
#   - DualSimplex
#   - DualSimplex + density filter
#   - CAM3
#   - Linseed
#   - TOAST
#
# Metric:
#   - RMSE (per cell type, normalized)
# ======================================================================

# -------------------------------
# Global settings
# -------------------------------
set.seed(3)

n_cell_types <- 15
n_genes      <- 10000
n_samples    <- 100
noise_sd     <- 5

dataset_name <- sprintf(
  "Simulated data (G=%d, N=%d, K=%d, sd=%g)",
  n_genes, n_samples, n_cell_types, noise_sd
)

# -------------------------------
# Libraries
# -------------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Deconvolution methods
library(DualSimplex)
library(CAM3)
library(TOAST)
library(linseed)

# Project utilities
source("R/simulation/create_simulation.R")
source("R/simulation/add_noise.R")

source("R/deconvolution/reference_free_methods.R")
source("R/filtering/iterative_density_filter.R")
source("R/metrics/metrics.R")
source("R/plotting/plotting_methods.R")

# ======================================================================
# Simulate data
# ======================================================================
message("Generating simulated dataset...")

sim <- create_simulation(
  n_genes = n_genes,
  n_samples = n_samples,
  n_cell_types = n_cell_types,
  with_marker_genes = FALSE
)

true_basis       <- sim$basis
true_proportions <- sim$proportions

sim <- sim %>% add_noise(noise_deviation = noise_sd)

expr <- sim$data

message("Input matrix dimensions: ",
        paste(dim(expr), collapse = " x "))

# ======================================================================
# Run reference-free deconvolution methods
# ======================================================================
results <- list()

# ----------------------------------------------------------------------
# 1. DualSimplex (no density filter)
# ----------------------------------------------------------------------
message("Running DualSimplex (no density filter)...")

res_ds <- deconvolve_dualsimplex(
  expr = expr,
  k = n_cell_types,
  use_density_filter = FALSE,
  verbose = TRUE
)

results[["DualSimplex"]] <- res_ds$proportions


# ----------------------------------------------------------------------
# 2. DualSimplex (with density filter)
# ----------------------------------------------------------------------
message("Running DualSimplex (with density filter)...")

res_ds_density <- deconvolve_dualsimplex(
  expr = expr,
  k = n_cell_types,
  use_density_filter = TRUE,
  density_threshold = 10,
  density_n_cell_types = n_cell_types,
  verbose = TRUE
)

results[["DualSimplex_density"]] <- res_ds_density$proportions


# ----------------------------------------------------------------------
# 3. CAM3
# ----------------------------------------------------------------------
message("Running CAM3...")

res_cam3 <- deconvolve_cam3(
  expr = expr,
  k = n_cell_types,
  thres.low = 0.3,
  thres.high = 0.95,
  radius.thres = 0.95,
  MG.num.thres = 20
)

results[["CAM3"]] <- res_cam3$proportions


# ----------------------------------------------------------------------
# 4. Linseed
# ----------------------------------------------------------------------
message("Running Linseed...")

res_linseed <- deconvolve_linseed(
  expr = expr,
  k = n_cell_types,
  topGenes = n_genes
)

results[["Linseed"]] <- res_linseed$proportions


# ----------------------------------------------------------------------
# 5. TOAST
# ----------------------------------------------------------------------
message("Running TOAST...")

res_toast <- deconvolve_toast(
  expr = expr,
  k = n_cell_types,
  nmarker = 4000
)

results[["TOAST"]] <- res_toast$proportions


# ======================================================================
# Metric computation
# ======================================================================
message("Computing RMSE metrics...")

estimated_list <- lapply(results, function(H) {
  H[1:nrow(true_proportions), , drop = FALSE]
})

metric_name <- "rmse_loss"

metric_df <- get_metric_values(
  estimated_list = estimated_list,
  true_proportions = true_proportions,
  metric = metric_name,
  per_row = TRUE,
  normalize = TRUE
)

# Order methods by median performance
method_order <- metric_df %>%
  group_by(method) %>%
  summarise(median_metric = median(metric_value)) %>%
  arrange(median_metric) %>%
  pull(method)

metric_df$method <- factor(metric_df$method, levels = method_order)

# ======================================================================
# Visualization
# ======================================================================
message("Generating comparison plot...")

colors <- brewer.pal(length(method_order), "Set1")
names(colors) <- method_order

comparisons <- list(
  c("DualSimplex_density", "DualSimplex"),
  c("DualSimplex_density", "CAM3"),
  c("DualSimplex_density", "Linseed"),
  c("DualSimplex_density", "TOAST")
)

p_rmse <- plot_single_boxplot_fixed(
  metric_results_for_datasets_and_methods = metric_df,
  comparisons = comparisons,
  title = dataset_name,
  reference_group = "DualSimplex",
  metric_column = "metric_value",
  metric_title = get_metric_plot_title(metric_name),
  test = "t.test",
  custom_colors = colors
)

print(p_rmse)

# ======================================================================
# Save results
# ======================================================================
output_dir <- file.path("results", "simulated")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(metric_df, file = file.path(output_dir, "rmse_metrics.rds"))
ggsave(
  filename = file.path(output_dir, "rmse_boxplot.pdf"),
  plot = p_rmse,
  width = 10,
  height = 6
)

message("Simulation benchmark completed successfully.")
