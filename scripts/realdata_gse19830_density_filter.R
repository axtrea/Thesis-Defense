#!/usr/bin/env Rscript
# ============================================================
# Density-based filtering and deconvolution on GSE19830
#
# Dataset: GSE19830 (rat tissue mixtures)
# Method: DualSimplex + iterative density filtering
# ============================================================


# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(biomaRt)
  library(rat2302.db)
  library(matrixStats)
  library(DualSimplex)
  library(linseed)
  library(ggplot2)
})


# ------------------------------------------------------------
# Global parameters
# ------------------------------------------------------------

n_ct <- 3
dataset_id <- "GSE19830"


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
  rat2302.db::rat2302.db,
  keys     = rownames(data_raw),
  columns  = c("SYMBOL", "GENETYPE"),
  keytype  = "PROBEID"
)

probe_mapping <- probe_mapping[
  probe_mapping$PROBEID %in% rownames(data_raw) &
    !is.na(probe_mapping$SYMBOL),
]

# One probe → one gene (dataset-specific cleanup)
probe_mapping <- probe_mapping[!duplicated(probe_mapping$PROBEID), ]

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

component_names <- c("Liver", "Brain", "Lung")
pdata <- pData(gse)

parsed_proportions <- strsplit(pdata$source_name_ch1, split = "/")
parsed_proportions <- lapply(parsed_proportions, function(x) {
  strtoi(gsub("\\D*(\\d+)\\D*", "\\1", x))
})

true_proportions <- matrix(
  unlist(parsed_proportions),
  ncol = n_ct,
  byrow = TRUE
) / 100

rownames(true_proportions) <- rownames(pdata)
colnames(true_proportions) <- component_names
true_proportions <- t(true_proportions)


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


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

message("Data dimensions:")
message("  Expression matrix: ", paste(dim(data_raw), collapse = " x "))
message("  True proportions:  ", paste(dim(true_proportions), collapse = " x "))
message("  True basis:        ", paste(dim(true_basis), collapse = " x "))


# ------------------------------------------------------------
# Initialize DualSimplex
# ------------------------------------------------------------

message("Initializing DualSimplex")

coding_gene_info <- probe_mapping[
  probe_mapping$GENETYPE == "protein-coding",
]

gene_anno <- list(
  CODING = coding_gene_info$SYMBOL
)

dso <- DualSimplexSolver$new()
dso$set_data(
  data_raw,
  gene_anno_lists = gene_anno,
  max_dim = 20
)

dso$basic_filter(
  remove_true_cols_default = c(),
  keep_true_cols = c("CODING")
)

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

set.seed(1)
dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()


# ------------------------------------------------------------
# Normalize estimated proportions
# ------------------------------------------------------------

cur_H <- solution$H

normalized_H <- apply(cur_H, 2, function(x) {
  x / sum(x)
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


message("Analysis finished.")
