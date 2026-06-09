# ============================================================
# Gaussian Noise Injection for Simulated Gene Expression Data
# ============================================================
# Author: Dmitrii Reshetin
#
# Description:
#   Adds Gaussian noise to simulated gene expression data.
#   Supports additive and multiplicative noise models with
#   optional protection of selected genes and samples.
# ============================================================


# ------------------------------------------------------------
# Add Gaussian noise to simulation object
# ------------------------------------------------------------
add_gaussian_noise <- function(
    simulation,
    noise_deviation,
    additive = TRUE,
    protect_genes = integer(0),
    protect_samples = integer(0)
) {
  
  if (is.null(simulation$data)) {
    stop("simulation$data is missing.")
  }
  
  data <- simulation$data
  
  # ----------------------------------------------------------
  # Generate Gaussian noise
  # ----------------------------------------------------------
  if (additive) {
    
    # Additive Gaussian noise (log2 scale)
    noise <- matrix(
      stats::rnorm(
        n = length(data),
        mean = 0,
        sd = noise_deviation
      ),
      nrow = nrow(data),
      ncol = ncol(data)
    )
    
    noise <- 2^noise
    
    # Protect selected genes and samples
    if (length(protect_genes) > 0) {
      noise[protect_genes, ] <- 0
    }
    if (length(protect_samples) > 0) {
      noise[, protect_samples] <- 0
    }
    
    noisy_data <- data + noise
    
  } else {
    
    # Multiplicative Gaussian noise (log2 scale)
    noise <- matrix(
      stats::rnorm(
        n = length(data),
        mean = 1,
        sd = noise_deviation
      ),
      nrow = nrow(data),
      ncol = ncol(data)
    )
    
    noise <- 2^noise
    
    # Protect selected genes and samples
    if (length(protect_genes) > 0) {
      noise[protect_genes, ] <- 1
    }
    if (length(protect_samples) > 0) {
      noise[, protect_samples] <- 1
    }
    
    noisy_data <- data * noise
  }
  
  # ----------------------------------------------------------
  # Ensure non-negativity
  # ----------------------------------------------------------
  noisy_data[noisy_data < 0] <- 0
  
  # ----------------------------------------------------------
  # Update simulation object
  # ----------------------------------------------------------
  simulation$data <- noisy_data
  simulation$gaussian_noise_params <- list(
    noise_deviation    = noise_deviation,
    additive           = additive,
    protected_genes    = protect_genes,
    protected_samples = protect_samples
  )
  
  return(simulation)
}
