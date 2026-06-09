#!/usr/bin/env Rscript
# ============================================================
# Density-based filtering and deconvolution on GSE11058
#
# Dataset: GSE11058 (cell line mixtures)
# Method: DualSimplex + iterative density filtering
# ============================================================


# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(biomaRt)
  library(hgu133plus2.db)
  library(matrixStats)
  library(DualSimplex)
  library(linseed)
  library(ggplot2)
})


# ------------------------------------------------------------
# Global parameters
# ------------------------------------------------------------

n_ct <- 4
dataset_id <- "GSE11058"


# ------------------------------------------------------------
# Load and preprocess GEO data
# ------------------------------------------------------------

message("Loading GEO dataset: ", dataset_id)

gse <- getGEO(dataset_id, AnnotGPL = TRUE)[[1]]
data_raw <- exprs(gse)

# Linearize expression values
data_raw <- linearize_dataset(data_raw)


# ------------------------------------------------------------
# Probe → gene mapping and gene-level aggregation
# ------------------------------------------------------------

message("Mapping probes to gene symbols")

probe_mapping <- biomaRt::select(
  hgu133plus2.db::hgu133plus2.db,
  keys     = rownames(data_raw),
  columns  = c("SYMBOL"),
  keytype  = "PROBEID"
)

probe_mapping <- probe_mapping[
  probe_mapping$PROBEID %in% rownames(data_raw) &
    !is.na(probe_mapping$SYMBOL),
]

data_raw <- data_raw[probe_mapping$PROBEID, ]
rownames(data_raw) <- probe_mapping$SYMBOL

# Collapse probes → genes using median
data_raw <- tapply(
  data_raw,
  list(
    row.names(data_raw)[row(data_raw)],
    colnames(data_raw)[col(data_raw)]
  ),
  FUN = median
)


# ------------------------------------------------------------
# Construct ground truth proportions (H)
# ------------------------------------------------------------

message("Constructing true proportions (H)")

pdata <- pData(gse)

component_names <- c("Jurkat", "IM-9", "Raji", "THP-1")

mix_a <- c(2.5, 1.25, 2.5, 3.75); mix_a <- mix_a / sum(mix_a)
mix_b <- c(0.5, 3.17, 4.75, 1.58); mix_b <- mix_b / sum(mix_b)
mix_c <- c(0.1, 4.95, 1.65, 3.3);  mix_c <- mix_c / sum(mix_c)
mix_d <- c(0.02, 3.33, 3.33, 3.33); mix_d <- mix_d / sum(mix_d)

proportions <- matrix(0, nrow = ncol(data_raw), ncol = n_ct)

proportions[1:3,  1] <- 1
proportions[4:6,  2] <- 1
proportions[7:9,  3] <- 1
proportions[10:12,4] <- 1

proportions[13:15,] <- mix_a
proportions[16:18,] <- mix_b
proportions[19:21,] <- mix_c
proportions[22:24,] <- mix_d

rownames(proportions) <- rownames(pdata)
colnames(proportions) <- component_names

true_proportions <- t(proportions)


# ------------------------------------------------------------
# Construct ground truth basis (W)
# ------------------------------------------------------------

message("Constructing true basis (W)")

components <- lapply(seq_len(n_ct), function(i) {
  cols <- colnames(true_proportions[, true_proportions[i, ] == 1])
  mat  <- data_raw[, cols]
  vec  <- as.matrix(rowMedians(mat))
  rownames(vec) <- rownames(mat)
  vec
})

true_basis <- do.call(cbind, components)
colnames(true_basis) <- component_names

true_marker_list <- get_signature_markers(
  true_basis,
  n_marker_genes = 100
)


# ------------------------------------------------------------
# Initialize DualSimplex
# ------------------------------------------------------------

message("Initializing DualSimplex")

dso <- DualSimplexSolver$new()
dso$set_data(data_raw)

# Remove non-coding genes (RPL, LOC, ORF, SNOR)
dso$basic_filter()

dso$project(n_ct)

message("Dimensions after basic filtering: ",
        paste(dim(dso$get_data()), collapse = " x "))


# ------------------------------------------------------------
# Density-based iterative filtering
# ------------------------------------------------------------

message("Running iterative density filter")

dso <- iterative_density_filter(
  dso,
  threshold     = 1,
  n_cell_types  = n_ct,
  genes         = TRUE
)

dso$project(n_ct)


# ------------------------------------------------------------
# Optimization
# ------------------------------------------------------------

message("Running optimization")

dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()


# ------------------------------------------------------------
# Normalize estimated proportions
# ------------------------------------------------------------

cur_H <- solution$H

normalized_H <- apply(cur_H, 2, function(x) {
  x <- x / sum(x)
  x
})


# ------------------------------------------------------------
# Evaluation
# ------------------------------------------------------------

message("Evaluating results")

ptp <- coerce_pred_true_props(normalized_H, true_proportions)
ptb <- coerce_pred_true_basis(
  solution$W,
  true_basis[rownames(solution$W), ]
)

plot_ptp_scatter(ptp)
plot_ptb_scatter(ptb)

linseed::plotProportions(
  as.data.frame(ptp[[1]]),
  as.data.frame(ptp[[2]]),
  pnames = c("predicted", "true"),
  point_size = 1,
  line_size  = 0.7
) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )


# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

dso$plot_negative_basis_change()
dso$plot_negative_proportions_change()

message("Analysis finished.")
