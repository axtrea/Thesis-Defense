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

# Функция для добавления отрицательного биномиального шума
add_negbinom_noise <- function(simulation, 
                               dispersion = 0.1,
                               protect_genes = c(), 
                               protect_samples = c()) {
  
  data <- simulation$data
  noisy_data <- matrix(0, nrow = nrow(data), ncol = ncol(data))
  
  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      if (!(i %in% protect_genes) && !(j %in% protect_samples)) {
        mu <- max(data[i, j], 0.1)  # Защита от нулей
        # size параметр связан с dispersion: var = mu + mu^2 * dispersion
        size <- 1 / dispersion
        noisy_data[i, j] <- stats::rnbinom(1, mu = mu, size = size)
      } else {
        noisy_data[i, j] <- data[i, j]
      }
    }
  }
  
  simulation$data <- noisy_data
  simulation$negbinom_noise_params <- list(
    dispersion = dispersion,
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
sim <- sim %>% add_negbinom_noise(dispersion = 1)  

data_raw <- sim$data
true_basis <- sim$basis
true_proportions <- sim$proportions

#Create dso object
dso <- DualSimplexSolver$new()
data_raw = getNonnegativeLowRankApproximationWithTangentMethod(data_raw, 3, 50, left=0)
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