# ============================================================
# kNN-based Density Annotation and Iterative Density Filter
# ============================================================
# Description:
#   Computes local density using kNN radius estimation and
#   applies an iterative density-based filtering procedure
#   in projected gene expression space.
# ============================================================

# ------------------------------------------------------------
# Density annotation
# ------------------------------------------------------------

add_knn_density_annotation <- function(
  eset,
  proj,
  genes = TRUE,
  k = 10,
  radius = NULL
) {

  anno <- get_anno(eset, genes)

  coords <- if (genes) proj$X else proj$Omega

  if (is.null(radius)) {
    knn <- dbscan::kNN(coords, k = k)
    radius <- median(knn$dist[, k], na.rm = TRUE)
  }

  nn <- dbscan::frNN(coords, eps = radius)

  density <- vapply(nn$id, length, integer(1))
  mean_dist <- vapply(nn$dist, mean, numeric(1))

  names(density) <- rownames(coords)
  names(mean_dist) <- rownames(coords)

  anno$density <- density[rownames(anno)]
  anno$mean_nn_distance <- mean_dist[rownames(anno)]

  set_anno(anno, eset, genes)
}


# ------------------------------------------------------------
# Iterative density filter
# ------------------------------------------------------------

iterative_density_filter <- function(
  dso,
  threshold = 1,
  n_cell_types,
  genes = FALSE,
  k = 10,
  max_iter = 100,
  verbose = TRUE
) {

  iter <- 1

  repeat {

    # 1. Add density annotation
    eset <- add_knn_density_annotation(
      eset  = dso$get_data(),
      proj  = dso$st$proj,
      genes = genes,
      k     = k
    )
    dso$set_data(eset)

    # 2. Filter by density
    filtered <- threshold_filter(
      eset       = eset,
      feature    = "density",
      threshold  = threshold,
      genes      = genes,
      keep_lower = FALSE
    )

    n_before <- if (genes) nrow(eset) else ncol(eset)
    n_after  <- if (genes) nrow(filtered) else ncol(filtered)
    removed  <- n_before - n_after

    if (verbose) {
      message(sprintf(
        "Iteration %d | removed %d objects (density < %.2f)",
        iter, removed, threshold
      ))
    }

    # 3. Stopping criteria
    if (removed == 0) break
    if (n_after == 0) {
      warning("All objects removed during density filtering.")
      break
    }
    if (iter >= max_iter) {
      warning("Maximum number of iterations reached.")
      break
    }

    # 4. Update and reproject
    dso$set_data(filtered)
    dso$project(n_cell_types)

    iter <- iter + 1
  }

  invisible(dso)
}
