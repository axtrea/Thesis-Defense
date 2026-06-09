#!/usr/bin/env Rscript
# ============================================================
# Density-based filtering and deconvolution on GSE64655
#
# Dataset: Hoek purified bulk RNA-seq
# Cell types: 6
# Method: DualSimplex + iterative density filtering
#
# ============================================================


# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(DualSimplex)
  library(Biobase)
  library(ggplot2)
})


# ------------------------------------------------------------
# Global parameters
# ------------------------------------------------------------

n_ct <- 6

tpm_file  <- "D:/denoising/validation/hoek_purified/hoek_purified_tpm.rds"
facs_file <- "D:/denoising/validation/hoek_purified/hoek_purified_facs.rds"


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

message("Loading Hoek purified dataset (GSE64655)")

bulk_expr <- readRDS(tpm_file)
fin_facs  <- readRDS(facs_file)


# ------------------------------------------------------------
# Expression preprocessing
# ------------------------------------------------------------

bulk_expr <- as.matrix(bulk_expr)
bulk_expr <- linearize_dataset(bulk_expr)
bulk_expr <- remove_zero_rows(bulk_expr)

input_matrix     <- bulk_expr
true_proportions <- as.matrix(fin_facs)


message("Expression matrix dimensions: ",
        paste(dim(input_matrix), collapse = " x "))
message("True proportions dimensions: ",
        paste(dim(true_proportions), collapse = " x "))


# ------------------------------------------------------------
# Initialize DualSimplex
# ------------------------------------------------------------

dso <- DualSimplexSolver$new()
dso$set_data(
  input_matrix,
  max_dim = 30,
  sinkhorn_tol = 1e-17
)

dso$project(n_ct)

dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)


# ------------------------------------------------------------
# Basic MAD-based prefiltering (as in Hoek benchmark)
# ------------------------------------------------------------

dummy_threshold <- 0.05

data <- dso$get_data()
anno <- get_anno(data)

anno$PASS_FILTER <- FALSE
anno[anno$log_mad > dummy_threshold, ]$PASS_FILTER <- TRUE

data <- set_anno(anno, data)

plot_feature(
  data,
  feature = "log_mad",
  col_by  = "PASS_FILTER"
)

dso$basic_filter(log_mad_gt = dummy_threshold)


# ------------------------------------------------------------
# Density annotation (before filtering)
# ------------------------------------------------------------

add_annov2(n_ct = n_ct, genes = TRUE)

density_values <- dso[["st"]][["data"]]@featureData@data[["density"]]
df <- data.frame(density = density_values)

ggplot(df, aes(x = density)) +
  geom_histogram(
    fill = "lightblue",
    color = "white",
    bins = 200
  ) +
  geom_vline(
    xintercept = 1,
    color = "red",
    linetype = "dashed",
    size = 1
  ) +
  labs(
    title = "Gene density distribution",
    x = "Density",
    y = "Frequency"
  ) +
  theme_minimal()


# ------------------------------------------------------------
# Iterative density filtering
# ------------------------------------------------------------

message("Running iterative density filter")

dso <- iterative_density_filter(
  dso,
  threshold    = 1,
  n_cell_types = n_ct,
  genes        = TRUE
)

dso$project(n_ct)

dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)


# ------------------------------------------------------------
# Optimization
# ------------------------------------------------------------

set.seed(1)

dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()


# ------------------------------------------------------------
# Evaluation
# ------------------------------------------------------------

ptp <- coerce_pred_true_props(solution$H, true_proportions)

plot_correlation_matrix(ptp[[1]], ptp[[2]])
plot_ptp_lines(ptp)
plot_ptp_scatter(ptp)

rmse <- rmse_loss_function(ptp[[1]], ptp[[2]])
message("RMSE: ", round(rmse, 4))


message("Analysis finished.")
