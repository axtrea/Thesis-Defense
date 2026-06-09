# ============================================================
# Uniform Noise Injection for Gene Expression Data
# ============================================================
# Description:
#   Adds uniform noise to simulated gene expression data.
#   Supports additive and multiplicative noise models.
# ============================================================


# ------------------------------------------------------------
# Add Uniform noise to simulation object
# ------------------------------------------------------------
add_uniform_noise <- function(
    simulation,
    noise_range = c(-0.5, 0.5),
    additive = TRUE,
    protect_genes = integer(0),
    protect_samples = integer(0)
) {
  
  if (is.null(simulation$data)) {
    stop("simulation$data is missing.")
  }
  
  if (length(noise_range) != 2 || noise_range[1] >= noise_range[2]) {
    stop("noise_range must be a numeric vector of length 2: c(min, max).")
  }
  
  data <- simulation$data
  
  # ----------------------------------------------------------
  # Generate Uniform noise
  # ----------------------------------------------------------
  if (additive) {
    
    # Additive uniform noise (log2 scale)
    noise <- matrix(
      stats::runif(
        n = length(data),
        min = noise_range[1],
        max = noise_range[2]
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
    
    # Multiplicative uniform noise (log2 scale)
    noise <- matrix(
      stats::runif(
        n = length(data),
        min = 1 + noise_range[1],
        max = 1 + noise_range[2]
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
  simulation$uniform_noise_params <- list(
    noise_range        = noise_range,
    type               = if (additive) "additive" else "multiplicative",
    protected_genes    = protect_genes,
    protected_samples = protect_samples
  )
  
  return(simulation)
}
