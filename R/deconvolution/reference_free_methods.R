# ---------------------------------------------------------------------
# Reference-free deconvolution methods
# ---------------------------------------------------------------------
# This file contains wrappers for popular reference-free deconvolution
# algorithms with a unified interface.
# ---------------------------------------------------------------------


#' DualSimplex deconvolution with optional iterative density filtering
#'
#' @param expr Numeric matrix (genes x samples)
#' @param k Integer, number of cell types
#' @param use_density_filter Logical, apply iterative density filter
#' @param density_threshold Numeric density cutoff
#' @param density_n_cell_types Integer, number of cell types for projection
#' @param max_density_iter Integer, max filtering iterations
#' @param max_sinkhorn_iter Integer, max Sinkhorn iterations
#' @param sinkhorn_tol Numeric tolerance for Sinkhorn
#' @param verbose Logical
#'
#' @return List with elements:
#' \item{object}{DualSimplexSolver object}
#' \item{proportions}{Estimated proportions matrix}
#'
#' @export
deconvolve_dualsimplex <- function(
  expr,
  k,
  use_density_filter = TRUE,
  density_threshold = 1,
  density_n_cell_types = NULL,
  max_density_iter = 100,
  max_sinkhorn_iter = 300,
  sinkhorn_tol = 1e-17,
  verbose = TRUE
) {

  if (is.null(density_n_cell_types)) {
    density_n_cell_types <- k
  }

  dso <- DualSimplexSolver$new()
  dso$set_data(
    expr,
    max_sinkhorn_iterations = max_sinkhorn_iter,
    max_dim = min(ncol(expr), 30),
    sinkhorn_tol = sinkhorn_tol
  )
  dso$project(k)

  if (verbose) {
    message("DualSimplex: initial data size = ",
            paste(dim(dso$st$data), collapse = " x "))
  }

  if (use_density_filter) {

    if (verbose) {
      message("Applying iterative density filter...")
    }

    dso$basic_filter(log_mad_gt = 0.05)

    dso <- iterative_density_filter(
      dso,
      threshold = density_threshold,
      n_cell_types = density_n_cell_types,
      genes = TRUE,
      max_iter = max_density_iter,
      verbose = verbose
    )

    if (verbose) {
      message("After filtering: data size = ",
              paste(dim(dso$st$data), collapse = " x "))
    }
  }

  dso$project(k)
  dso$init_solution("random")
  dso$default_optimization()

  sol <- dso$finalize_solution()

  list(
    object = dso,
    proportions = sol$H
  )
}


#' TOAST reference-free deconvolution
#'
#' @param expr Numeric matrix (genes x samples)
#' @param k Integer, number of cell types
#' @param nmarker Integer, number of marker genes
#'
#' @return List with elements object and proportions
#' @export
deconvolve_toast <- function(expr, k, nmarker = 4000) {

  data_raw <- as.matrix(expr)

  refinx <- findRefinx(data_raw, nmarker = nmarker)
  Y <- data_raw[refinx, , drop = FALSE]
  Y <- sweep(Y, 2, colSums(Y), "/")

  out <- myRefFreeCellMix(
    Y,
    mu0 = myRefFreeCellMixInitialize(Y, K = k),
    iters = 500,
    verbose = 0
  )

  proportions <- t(out$Omega)
  colnames(proportions) <- colnames(data_raw)
  rownames(proportions) <- paste0("toast_cell_type_", seq_len(k))

  list(
    object = out,
    proportions = proportions
  )
}


#' CAM3 reference-free deconvolution
#'
#' @param expr Numeric matrix (genes x samples)
#' @param k Integer, number of cell types
#'
#' @return List with elements object and proportions
#' @export
deconvolve_cam3 <- function(
  expr,
  k,
  thres.low = 0.3,
  thres.high = 0.95,
  radius.thres = 0.95,
  MG.num.thres = 20
) {

  data_raw <- as.matrix(expr)
  Y <- sweep(data_raw, 2, colSums(data_raw), "/")

  cam <- CAM3Run(
    Y,
    K = k,
    dim.rdc = k,
    cluster.num = 120,
    thres.low = thres.low,
    thres.high = thres.high,
    radius.thres = radius.thres,
    MG.num.thres = MG.num.thres
  )

  proportions <- t(cam@ASestResult[[1]]@Aest)
  colnames(proportions) <- colnames(data_raw)
  rownames(proportions) <- paste0("cam3_cell_type_", seq_len(k))

  list(
    object = cam,
    proportions = proportions
  )
}


#' Linseed reference-free deconvolution
#'
#' @param expr Numeric matrix (genes x samples)
#' @param k Integer, number of cell types
#' @param topGenes Integer
#'
#' @return List with elements object and proportions
#' @export
deconvolve_linseed <- function(expr, k, topGenes = 8000) {

  data_raw <- as.matrix(expr)

  lo <- LinseedObject$new(data_raw, topGenes = topGenes)
  lo$calculatePairwiseLinearity()
  lo$calculateSpearmanCorrelation()
  lo$calculateSignificanceLevel(100)
  lo$filterDatasetByPval(0.01)
  lo$setCellTypeNumber(k)
  lo$project("filtered")
  lo$smartSearchCorners(dataset = "filtered", error = "norm")

  proportions <- lo$proportions
  colnames(proportions) <- colnames(data_raw)
  rownames(proportions) <- paste0("linseed_cell_type_", seq_len(k))

  list(
    object = lo,
    proportions = proportions
  )
}


#' CellDistinguisher reference-free deconvolution
#'
#' @param expr Numeric matrix (genes x samples)
#' @param k Integer, number of cell types
#'
#' @return List with elements object and proportions
#' @export
deconvolve_celldistinguisher <- function(expr, k) {

  data_raw <- as.matrix(expr)
  expr_linear <- data_raw[rowSums(data_raw == 0) != ncol(data_raw), ]

  dist <- gecd_CellDistinguisher(
    expr_linear,
    genesymb = NULL,
    numCellClasses = k,
    minDistinguisherAlternatives = 20,
    maxDistinguisherAlternatives = 100,
    minAlternativesLengthsNormalized = 0.5,
    expressionQuantileForFilter = 0.995,
    expressionConcentrationRatio = 0.333,
    verbose = 0
  )

  dec <- gecd_DeconvolutionCellMix(
    expr_linear,
    dist$bestDistinguishers,
    method = "ssKL",
    maxIter = 5
  )

  proportions <- dec$sampleCompositions
  colnames(proportions) <- colnames(data_raw)
  rownames(proportions) <- paste0("cd_cell_type_", seq_len(k))

  list(
    object = dec,
    proportions = proportions
  )
}
