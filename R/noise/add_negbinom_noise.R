# ============================================================
# Negative Binomial Noise Injection for Gene Expression Data
# ============================================================
# Description:
#   Adds Negative Binomial (NB) noise to simulated gene
#   expression data. This noise model reflects overdispersed
#   count variability commonly observed in RNA-seq data.
#
#   The dispersion parameter controls variance according to:
#     Var(X) = mu + mu^2 * dispersion
#
#   Optional protection of genes and samples allows controlled
#   perturbation of the expression matrix.
# ============================================================


# ------------------------------------------------------------
# Add Negative Binomial noise to simulation object
# ------------------------------------------------------------
add_negbinom_noise <- function(
    simulation,
    dispersion = 0.1,
    protect_genes = integer(0),
    protect_samples = integer(0),
    min_mu = 0.1
) {
  
  if (is.null(simulation$data)) {
    stop("simulation$data is missing.")
  }
  
  data <- simulation$data
  noisy_data <- matrix(
    0,
    nrow = nrow(data),
    ncol = ncol(data),
    dimnames = dimnames(data)
  )
  
  # NB size parameter from dispersion
  size <- 1 / dispersion
  
  # ----------------------------------------------------------
  # Noise injection
  # ----------------------------------------------------------
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(ncol(data))) {
      
      if (i %in% protect_genes || j %in% protect_samples) {
        noisy_data[i, j] <- data[i, j]
        next
      }
      
      # Mean expression (protected from zeros)
      mu <- max(data[i, j], min_mu)
      
      noisy_data[i, j] <- stats::rnbinom(
        n = 1,
        mu = mu,
        size = size
      )
    }
  }
  
  # ----------------------------------------------------------
  # Update simulation object
  # ----------------------------------------------------------
  simulation$data <- noisy_data
  simulation$negbinom_noise_params <- list(
    dispersion         = dispersion,
    size_parameter     = size,
    min_mu             = min_mu,
    protected_genes    = protect_genes,
    protected_samples = protect_samples
  )
  
  return(simulation)
}
