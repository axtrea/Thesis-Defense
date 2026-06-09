# ============================================================
# Poisson Noise Injection for Gene Expression Data
# ============================================================
# Description:
#   Adds Poisson noise to simulated gene expression data.
#   This model approximates count-based sampling noise
#   commonly observed in transcriptomics.
# ============================================================


# ------------------------------------------------------------
# Add Poisson noise to simulation object
# ------------------------------------------------------------
add_poisson_noise <- function(
    simulation,
    lambda_noise = 1,
    protect_genes = integer(0),
    protect_samples = integer(0),
    eps = 1e-8
) {
  
  if (is.null(simulation$data)) {
    stop("simulation$data is missing.")
  }
  
  data <- simulation$data
  
  # ----------------------------------------------------------
  # Numerical safety
  # ----------------------------------------------------------
  data_safe <- pmax(data, eps)
  
  # ----------------------------------------------------------
  # Generate Poisson noise
  # ----------------------------------------------------------
  poisson_noise <- matrix(
    stats::rpois(
      n = length(data_safe),
      lambda = lambda_noise
    ),
    nrow = nrow(data_safe),
    ncol = ncol(data_safe)
  )
  
  poisson_noise <- 2^poisson_noise
  
  # ----------------------------------------------------------
  # Protect selected genes and samples
  # ----------------------------------------------------------
  if (length(protect_genes) > 0) {
    poisson_noise[protect_genes, ] <- 0
  }
  if (length(protect_samples) > 0) {
    poisson_noise[, protect_samples] <- 0
  }
  
  # ----------------------------------------------------------
  # Apply noise
  # ----------------------------------------------------------
  noisy_data <- data + poisson_noise
  
  # Ensure non-negativity
  noisy_data[noisy_data < 0] <- 0
  
  # ----------------------------------------------------------
  # Update simulation object
  # ----------------------------------------------------------
  simulation$data <- noisy_data
  simulation$poisson_noise_params <- list(
    lambda_noise     = lambda_noise,
    eps              = eps,
    protected_genes  = protect_genes,
    protected_samples = protect_samples
  )
  
  return(simulation)
}
