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

# Function to calculate distance from a point to a triangle
distance_to_triangle <- function(point, triangle_vertices) {
  distances <- sqrt(rowSums((t(t(triangle_vertices) - point))^2))
  return(min(distances))
}

find_k_nearest_genes <- function(projected_points, true_vertices, k = 500) {
  # Calculate distances for all points
  distances <- apply(projected_points, 1, function(point) {
    distance_to_triangle(point, true_vertices)
  })
  
  nearest_indices <- order(distances)[1:k]
  
  return(list(
    indices = nearest_indices,
    distances = distances[nearest_indices],
    all_distances = distances
  ))
}

# util method to add triangles
add_points_to_plot <- function(plt, points, with_lines=T,  with_points =T,
                               color_points="purple", color_lines="blue", line_type="dashed",
                               points_label = "True solution") {
  points <- as.data.frame(points)
  points$color <- points_label
  
  if (with_lines) {
    plt <-  plt + 
      geom_polygon(
        data = points,
        size = 1,
        fill = NA,
        color = color_lines,
        linetype = line_type,
        aes(fill = points_label)
      ) 
  }
  if (with_points) {
    plt <- plt +
      geom_point(
        data = points,
        color = color_points,
        size = 3,
        aes(fill = points_label)
      )  
  }
  return(plt)
}

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

set.seed(3)
sim <- create_simulation(n_genes = 10000,
                         n_samples = 100,
                         n_cell_types = 3,
                         with_marker_genes = FALSE)

data_raw <- sim$data
true_basis <- sim$basis
true_proportions <- sim$proportions

dso <- DualSimplexSolver$new()
dso$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
dso$project(3)

cropped_basis <-  true_basis[rownames(dso$st$data),]

extended_scaling_result <- extended_sinkhorn_scale(V = cropped_basis %*% true_proportions,
                                                   W=cropped_basis[rownames(dso$st$data),],
                                                   H=true_proportions[, colnames(dso$st$data)],
                                                   n_iter = dso$st$scaling$iterations)

# true points
true_X_points <-   (extended_scaling_result$W_row %*% extended_scaling_result$H_row) %*% t(dso$st$proj$meta$R)
true_solution_coordinates <-  dso$get_coordinates_from_external_matrices(true_basis, true_proportions)

# ========== FIND 500 NEAREST GENES ==========
clean_projected_points <- dso$st$proj$X[, 2:3]

# Find 500 nearest genes to the true solution triangle
nearest_genes <- find_k_nearest_genes(
  projected_points = clean_projected_points,
  true_vertices = true_solution_coordinates$X[, 2:3],
  k = 500
)

# Store the indices of nearest genes
nearest_gene_indices <- nearest_genes$indices

# Subset the true_X_points to only include these 500 genes
true_X_points_filtered <- true_X_points[nearest_gene_indices, ]

# ========== ADD NOISE AND PROJECT NOISY DATA ==========

# 1. Add noise to the original simulation
sim_noisy <- sim %>% add_gaussian_noise(noise_deviation = 0.5, additive = TRUE)

# 2. Create new DualSimplexSolver with noisy data
dso_noisy <- DualSimplexSolver$new()
dso_noisy$set_data(sim_noisy$data, max_dim = 30, sinkhorn_tol = 1e-17)
dso_noisy$project(3)

# 3. Get noisy projection points
projection_plot_noisy <- plot_projection_points(dso_noisy$st$proj, use_dims = (2:3), spaces = "X")
projection_plot_noisy$data <- projection_plot_noisy$data[nearest_gene_indices, ]

original_plot_X_noisy <- projection_plot_noisy %>% add_points_to_plot(true_solution_coordinates$X,  
                                                                      color_lines ="blue",
                                                                      color_points = "red", points_label = NULL)

connector_data_X_true <- data.frame(
  x_start = dso$st$proj$X[,2], 
  y_start = dso$st$proj$X[,3],
  x_end   = true_X_points[,2], 
  y_end   =  true_X_points[,3]
)

original_plot_X_noisy + geom_segment(data = connector_data_X_true, 
                                     aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                                     color = "blue", linetype = "dashed")+ geom_point(data=true_X_points_filtered, color= "green")
