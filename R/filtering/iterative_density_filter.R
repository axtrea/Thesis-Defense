# ============================================================
# Iterative Density Filter (DualSimplex-compatible)
# ============================================================
# Author: Dmitrii Reshetin
#
# Description:
#   Iterative density-based filtering for gene expression data
#   integrated with DualSimplexSolver. Density is computed
#   using a fixed MAD-based radius in projected space and
#   iteratively removes low-density outliers.
# ============================================================

library(dbscan)

# ------------------------------------------------------------
# Density annotation (MAD-based fixed radius)
# ------------------------------------------------------------
add_density_annotation <- function(
    eset,
    proj,
    genes = TRUE,
    n_cell_types = 3,
    radius = NULL
) {
  
  anno <- get_anno(eset, genes)
  
  # Select projected coordinates
  coords <- if (genes) proj$X else proj$Omega
  
  # Robust radius estimation (MAD-based)
  if (is.null(radius)) {
    if (is.matrix(coords) || is.data.frame(coords)) {
      coords_mad <- if (ncol(coords) > 1) {
        coords[, 2:ncol(coords), drop = FALSE]
      } else {
        coords
      }
      radius <- stats::mad(coords_mad) * sqrt(n_cell_types)
    } else {
      radius <- stats::mad(coords) * sqrt(n_cell_types)
    }
  }
  
  # Fixed-radius nearest neighbors
  nn <- dbscan::frNN(coords, eps = radius)
  
  density <- vapply(nn$id, length, numeric(1))
  mean_dist <- vapply(nn$dist, function(x) mean(x, na.rm = TRUE), numeric(1))
  
  names(density) <- rownames(coords)
  names(mean_dist) <- rownames(coords)
  
  anno$density <- density[rownames(anno)]
  anno$mean_nn_distance <- mean_dist[rownames(anno)]
  
  eset <- set_anno(anno, eset, genes)
  return(eset)
}

# ------------------------------------------------------------
# Iterative Density Filter
# ------------------------------------------------------------
iterative_density_filter <- function(
    dso,
    threshold = 1,
    n_cell_types = 3,
    genes = TRUE,
    max_iter = 100,
    verbose = TRUE
) {
  
  iter <- 1
  
  repeat {
    
    if (verbose) {
      message(
        sprintf(
          "Iteration %d | computing density annotation",
          iter
        )
      )
    }
    
    # 1. Update density annotation
    eset <- add_density_annotation(
      eset = dso$get_data(),
      proj = dso$st$proj,
      genes = genes,
      n_cell_types = n_cell_types
    )
    
    # 2. Apply threshold filter
    filtered_eset <- threshold_filter(
      eset = eset,
      feature = "density",
      threshold = threshold,
      genes = genes,
      keep_lower = FALSE
    )
    
    # 3. Track changes
    n_before <- if (genes) nrow(eset) else ncol(eset)
    n_after  <- if (genes) nrow(filtered_eset) else ncol(filtered_eset)
    removed <- n_before - n_after
    
    if (verbose) {
      message(
        sprintf(
          "Iteration %d | removed %d objects (density < %g)",
          iter, removed, threshold
        )
      )
    }
    
    # --- stopping criteria ---
    if (removed == 0) {
      if (verbose) message("Convergence reached: no objects removed.")
      break
    }
    
    if (n_after == 0) {
      warning("All objects removed — stopping.")
      dso$set_data(filtered_eset)
      break
    }
    
    if (iter >= max_iter) {
      warning("Maximum iterations reached — stopping.")
      dso$set_data(filtered_eset)
      break
    }
    
    # 4. Update solver state
    dso$set_data(filtered_eset)
    dso$project(n_cell_types)
    
    iter <- iter + 1
  }
  
  invisible(dso)
}
