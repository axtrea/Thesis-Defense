library("DualSimplex")
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(svglite)
library(ggpubr)
library(ComplexHeatmap)
library(progress)
library(reshape)
library(plotly)
library(ggrastr)
library(linseed)
library(matrixStats)

# Функция для добавления гауссова шума (без изменений)
add_gaussian_noise <- function(simulation, noise_deviation, 
                               additive = TRUE, 
                               protect_genes = c(), 
                               protect_samples = c()) {
  
  data <- simulation$data
  
  if (additive) {
    # Аддитивный гауссов шум
    noise <- matrix(
      stats::rnorm(n = length(data), mean = 0, sd = noise_deviation),
      nrow = nrow(data),
      ncol = ncol(data)
    )
    noise <- 2^noise
    # Защита генов и образцов
    if (length(protect_genes) > 0) {
      noise[protect_genes, ] <- 0
    }
    if (length(protect_samples) > 0) {
      noise[, protect_samples] <- 0
    }
    
    noisy_data <- data + noise
    
  } else {
    # Мультипликативный гауссов шум
    noise <- matrix(
      stats::rnorm(n = length(data), mean = 1, sd = noise_deviation),
      nrow = nrow(data),
      ncol = ncol(data)
    )
    noise <- 2^noise
    # Защита генов и образцов
    if (length(protect_genes) > 0) {
      noise[protect_genes, ] <- 1
    }
    if (length(protect_samples) > 0) {
      noise[, protect_samples] <- 1
    }
    
    noisy_data <- data * noise
  }
  
  # Обеспечение неотрицательности
  noisy_data[noisy_data < 0] <- 0
  
  simulation$data <- noisy_data
  simulation$gaussian_noise_params <- list(
    noise_deviation = noise_deviation,
    additive = additive,
    protected_genes = protect_genes,
    protected_samples = protect_samples
  )
  
  simulation
}

#Create simulation
set.seed(3)
sim <- create_simulation(n_genes = 10000,
                         n_samples = 100,
                         n_cell_types = 3,
                         with_marker_genes = FALSE)
sim <- sim %>% add_gaussian_noise(noise_deviation = 5, additive = TRUE)

data_raw <- sim$data
true_basis <- sim$basis
true_proportions <- sim$proportions
data_raw = getNonnegativeLowRankApproximationWithTangentMethod(data_raw, 3, 50, left=0)
#Create dso object
dso <- DualSimplexSolver$new()
dso$set_data(data_raw$newX, max_dim = 30, sinkhorn_tol = 1e-17)
dso$project(3)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

#Solution
dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()
names(solution)
solution <- dso$get_solution()

#Plot
ptp <- coerce_pred_true_props(solution$H, true_proportions)
plot_ptp_scatter(ptp)
plot_ptp_lines(ptp)