library("DualSimplex")
library(dplyr)
library(reshape2) 
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(RColorBrewer)
library(dbscan)

compute_extremeness <- function(coords, k = 10) {
  # Центроид
  centroid <- colMeans(coords)
  
  # Расстояние до центроида
  dist_to_centroid <- sqrt(rowSums((coords - matrix(centroid, 
                                                    nrow(coords), 
                                                    ncol(coords), 
                                                    byrow = TRUE))^2))
  
  # kNN расстояния
  knn <- dbscan::kNN(coords, k = k)
  knn_dist <- rowMeans(knn$dist, na.rm = TRUE)
  
  # Нормализация
  ext1 <- scale(dist_to_centroid)[,1]
  ext2 <- scale(knn_dist)[,1]
  
  extremeness <- ext1 + ext2
  
  names(extremeness) <- rownames(coords)
  return(extremeness)
}

add_extremeness_anno <- function(dso, k = 10, genes = TRUE) {
  
  eset <- dso$get_data()
  
  # координаты
  coords <- if (genes) dso$st$proj$X else dso$st$proj$Omega
  
  extremeness <- compute_extremeness(coords, k = k)
  
  anno <- get_anno(eset, genes)
  anno$extremeness <- extremeness[rownames(anno)]
  
  eset <- set_anno(anno, eset, genes)
  
  dso$set_data(eset)
  
  return(dso)
}

compute_extremeness <- function(coords, k = 10) {
  # Центроид
  centroid <- colMeans(coords)
  
  # Расстояние до центроида
  dist_to_centroid <- sqrt(rowSums((coords - matrix(centroid, 
                                                    nrow(coords), 
                                                    ncol(coords), 
                                                    byrow = TRUE))^2))
  
  # kNN расстояния
  knn <- dbscan::kNN(coords, k = k)
  knn_dist <- rowMeans(knn$dist, na.rm = TRUE)
  
  # Нормализация
  ext1 <- scale(dist_to_centroid)[,1]
  ext2 <- scale(knn_dist)[,1]
  
  extremeness <- ext1 + ext2
  
  names(extremeness) <- rownames(coords)
  return(extremeness)
}

set.seed(3)
n_ct = 8
sim <- create_simulation(n_genes = 10000,
                         n_samples = 100,
                         n_cell_types = n_ct,
                         with_marker_genes = FALSE)

true_basis <- sim$basis
true_proportions <- sim$proportions


sim <- sim %>% add_noise(noise_deviation = 3)
data_raw <- sim$data
dso <- DualSimplexSolver$new()
dso$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
dso$project(n_ct)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

add_extremeness_anno(dso = dso, k = 10, genes = TRUE)

eset <- dso$get_data()   

# Извлекаем вектор extremeness
extremeness_values <- eset@featureData@data[["extremeness"]]

# Строим гистограмму с помощью ggplot2
library(ggplot2)

ggplot(data.frame(extremeness = extremeness_values), aes(x = extremeness)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Extremeness (genes)",
       x = "Extremeness score",
       y = "Frequency") +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 0, y = Inf, label = "zero", vjust = 2, color = "red")

eset <- dso$get_data()
data <- threshold_filter(eset, feature = "extremeness", threshold = 0, keep_lower = F)
dso$set_data(data)
dso$project(n_ct)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()
ptp <- coerce_pred_true_props(solution$H, true_proportions)
plot_ptp_scatter(ptp)
plot_ptp_lines(ptp)