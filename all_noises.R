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

#' Добавить пакетные эффекты (batch effects)
#'
#' @param simulation оригинальный объект симуляции
#' @param n_batches количество батчей
#' @param batch_strength сила пакетного эффекта (0-1)
#' @param multiplicative если TRUE - мультипликативный эффект, FALSE - аддитивный
#' @param protect_genes не модифицировать эти строки
#' @param protect_samples не модифицировать эти столбцы
#' @param seed seed для воспроизводимости
#' @return модифицированный объект симуляции
#' @export
add_batch_effects <- function(simulation, n_batches = 3,
                              batch_strength = 0.3,
                              multiplicative = TRUE,
                              protect_genes = c(), 
                              protect_samples = c(),
                              seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  data <- simulation$data
  n_samples <- ncol(data)
  
  # Назначаем образцы случайным батчам
  batch_assignments <- sample(1:n_batches, n_samples, replace = TRUE)
  
  # Генерируем эффекты для каждого батча
  if (multiplicative) {
    # Мультипликативные эффекты
    batch_factors <- stats::rnorm(n_batches, mean = 1, sd = batch_strength)
    batch_factors[batch_factors < 0] <- 0.1  # Защита от отрицательных
    
    noisy_data <- data
    for (b in 1:n_batches) {
      batch_samples <- which(batch_assignments == b)
      if (length(batch_samples) > 0) {
        noisy_data[, batch_samples] <- data[, batch_samples] * batch_factors[b]
      }
    }
  } else {
    # Аддитивные эффекты
    batch_offsets <- stats::rnorm(n_batches, mean = 0, sd = batch_strength * mean(data))
    
    noisy_data <- data
    for (b in 1:n_batches) {
      batch_samples <- which(batch_assignments == b)
      if (length(batch_samples) > 0) {
        noisy_data[, batch_samples] <- data[, batch_samples] + batch_offsets[b]
      }
    }
  }
  
  # Обеспечение неотрицательности
  noisy_data[noisy_data < 0] <- 0
  
  # Защита генов и образцов
  if (length(protect_genes) > 0) {
    noisy_data[protect_genes, ] <- data[protect_genes, ]
  }
  if (length(protect_samples) > 0) {
    noisy_data[, protect_samples] <- data[, protect_samples]
  }
  
  simulation$data <- noisy_data
  simulation$batch_effects_params <- list(
    n_batches = n_batches,
    batch_strength = batch_strength,
    multiplicative = multiplicative,
    batch_assignments = batch_assignments
  )
  simulation$protected_genes <- protect_genes
  simulation$protected_samples <- protect_samples
  
  simulation
}

#' Добавить Dropout шум (имитация нулей в scRNA-seq)
#'
#' @param simulation оригинальный объект симуляции
#' @param dropout_rate базовый уровень dropout (0-1)
#' @param mean_based если TRUE, вероятность dropout зависит от среднего выражения
#' @param protect_genes не модифицировать эти строки
#' @param protect_samples не модифицировать эти столбцы
#' @param seed seed для воспроизводимости
#' @return модифицированный объект симуляции
#' @export
add_dropout_noise <- function(simulation, dropout_rate = 0.1,
                              mean_based = TRUE,
                              protect_genes = c(), 
                              protect_samples = c(),
                              seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  data <- simulation$data
  
  if (mean_based) {
    # Вероятность dropout обратно пропорциональна уровню выражения
    # Модель: p(dropout) = 1 / (1 + exp(-(α - β*log(expression))))
    
    # Нормализуем данные
    log_data <- log1p(data)
    scaled_data <- (log_data - mean(log_data)) / sd(log_data)
    
    # Вероятность dropout
    dropout_probs <- 1 / (1 + exp(-(1 - 2 * scaled_data))) * dropout_rate
    
  } else {
    # Равномерная вероятность
    dropout_probs <- matrix(dropout_rate, 
                            nrow = nrow(data), 
                            ncol = ncol(data))
  }
  
  # Создаем маску dropout
  dropout_mask <- matrix(
    stats::runif(n = length(data)) > dropout_probs,
    nrow = nrow(data),
    ncol = ncol(data)
  )
  
  noisy_data <- data * dropout_mask
  
  # Защита генов и образцов
  if (length(protect_genes) > 0) {
    noisy_data[protect_genes, ] <- data[protect_genes, ]
  }
  if (length(protect_samples) > 0) {
    noisy_data[, protect_samples] <- data[, protect_samples]
  }
  
  simulation$data <- noisy_data
  simulation$dropout_noise_params <- list(
    dropout_rate = dropout_rate,
    mean_based = mean_based,
    zeros_proportion = sum(noisy_data == 0) / length(noisy_data)
  )
  simulation$protected_genes <- protect_genes
  simulation$protected_samples <- protect_samples
  
  simulation
}

#' Добавить Negative Binomial шум (для RNA-seq счетов)
#'
#' @param simulation оригинальный объект симуляции
#' @param dispersion параметр дисперсии (φ), чем больше - тем выше дисперсия
#' @param scaling_factor фактор масштабирования λ
#' @param protect_genes не модифицировать эти строки
#' @param protect_samples не модифицировать эти столбцы
#' @param seed seed для воспроизводимости
#' @return модифицированный объект симуляции
#' @export
add_negative_binomial_noise <- function(simulation, dispersion = 0.1,
                                        scaling_factor = 1,
                                        protect_genes = c(), 
                                        protect_samples = c(),
                                        seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  data <- simulation$data
  
  # Negative Binomial: Var = μ + φμ²
  # size = 1/φ
  size_param <- 1 / dispersion
  
  # Генерируем шум
  noisy_data <- matrix(
    stats::rnbinom(n = length(data),
                   mu = data * scaling_factor,
                   size = size_param),
    nrow = nrow(data),
    ncol = ncol(data)
  )
  
  # Защита генов и образцов
  if (length(protect_genes) > 0) {
    noisy_data[protect_genes, ] <- data[protect_genes, ]
  }
  if (length(protect_samples) > 0) {
    noisy_data[, protect_samples] <- data[, protect_samples]
  }
  
  simulation$data <- noisy_data
  simulation$nb_noise_params <- list(
    dispersion = dispersion,
    scaling_factor = scaling_factor,
    size = size_param
  )
  simulation$protected_genes <- protect_genes
  simulation$protected_samples <- protect_samples
  
  simulation
}

#' Добавить комбинированный биологический шум (рекомендуется)
#'
#' @param simulation оригинальный объект симуляции
#' @param dropout_rate уровень dropout
#' @param dispersion дисперсия Negative Binomial
#' @param batch_strength сила пакетных эффектов
#' @param protect_genes не модифицировать эти строки
#' @param protect_samples не модифицировать эти столбцы
#' @param seed seed для воспроизводимости
#' @return модифицированный объект симуляции
#' @export
add_combined_biological_noise <- function(simulation,
                                          dropout_rate = 0.05,
                                          dispersion = 0.2,
                                          batch_strength = 0.2,
                                          protect_genes = c(), 
                                          protect_samples = c(),
                                          seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Клонируем симуляцию для последовательного добавления шума
  sim_noisy <- simulation
  
  # 1. Пакетные эффекты (сначала, как в реальном эксперименте)
  sim_noisy <- add_batch_effects(sim_noisy,
                                 n_batches = 2,
                                 batch_strength = batch_strength,
                                 multiplicative = TRUE,
                                 seed = seed)
  
  # 2. Negative Binomial шум (техническая вариация)
  sim_noisy <- add_negative_binomial_noise(sim_noisy,
                                           dispersion = dispersion,
                                           seed = if(!is.null(seed)) seed + 1 else NULL)
  
  # 3. Dropout эффекты
  sim_noisy <- add_dropout_noise(sim_noisy,
                                 dropout_rate = dropout_rate,
                                 mean_based = TRUE,
                                 seed = if(!is.null(seed)) seed + 2 else NULL)
  
  # Защита генов и образцов (применяем в конце)
  if (length(protect_genes) > 0) {
    sim_noisy$data[protect_genes, ] <- simulation$data[protect_genes, ]
  }
  if (length(protect_samples) > 0) {
    sim_noisy$data[, protect_samples] <- simulation$data[, protect_samples]
  }
  
  sim_noisy$combined_noise_params <- list(
    dropout_rate = dropout_rate,
    dispersion = dispersion,
    batch_strength = batch_strength,
    steps = c("batch_effects", "negative_binomial", "dropout")
  )
  sim_noisy$protected_genes <- protect_genes
  sim_noisy$protected_samples <- protect_samples
  
  sim_noisy
}

#Simulation
set.seed(5)
sim_initial <- create_simulation(n_genes = 10000, n_samples = 100, n_cell_types = 3)
sim_combined <- sim_initial %>% 
  add_combined_biological_noise(dropout_rate = 0.5, dispersion = 10, seed = 5)

dso <- DualSimplexSolver$new()
dso$set_data(sim_combined$data, max_dim = 30, sinkhorn_tol = 1e-17)
dso$project(3)
dso$plot_svd(dims=1:15)
dso$plot_projected("zero_distance", "zero_distance", use_dims = list(2:3))

#LowRankApproximation
data_new = getNonnegativeLowRankApproximationWithTangentMethod(sim_combined$data, 3, 50, left=0)
dso$set_data(data_new$newX, max_dim = 30, sinkhorn_tol = 1e-17)
dso$project(3)
dso$plot_svd(dims=1:15)
dso$plot_projected("zero_distance", "zero_distance", use_dims = list(2:3))

#MAD threshold
dummy_threshold <- 0.1
data <- dso$get_data()
anno <- get_anno(data)
anno$PASS_FILTER <- FALSE
anno[anno$log_mad > dummy_threshold,]$PASS_FILTER <- TRUE
data <- set_anno(anno, data)
plot_feature(data, feature = "log_mad", col_by ='PASS_FILTER')

dso$basic_filter(log_mad_gt = dummy_threshold)
dso$project(3)

svd_plot_2 <- dso$plot_svd(1:10) + theme_minimal(base_size = 8)
points_2 <- dso$plot_projected("black", "black",use_dims = 2:3, show_plots=F) +  theme_minimal(base_size = 8)
plotlist = list(points_2, svd_plot_2)
cowplot::plot_grid(plotlist=plotlist, rel_widths = c(0.66, 0.33))

#Denoising
distance_dist_plot_noizy <-
  plot_feature(dso$get_data(), feature = "plane_distance") + theme_minimal(base_size = 15)
color = fData(dso$get_data())$plane_distance
X_points_noizy <-
  plot_projection_points(
    dso$st$proj,
    use_dims = 2:3,
    color = color,
    color_name = "plane_distance",
    spaces = c("X"),
    pt_size = 2
  ) + theme_minimal(base_size = 15)

distances_distribution <-
  plot_feature_pair(dso$get_data(), "plane_distance", "zero_distance", T, size = 0.1)


plotlist = list(X_points_noizy,
                distance_dist_plot_noizy,
                distances_distribution)

cowplot::plot_grid(
  plotlist = plotlist,
  nrow = 1,
  rel_widths = c(0.4, 0.3, 0.3)
)


dso$distance_filter(plane_d_lt = 0.5, zero_d_lt = 0.15, genes = T)
dso$project(3)

current_svd_plot <- dso$plot_svd(1:10) + theme_minimal(base_size = 8)
current_points_plot <-
  dso$plot_projected("plane_distance",
                     "plane_distance",
                     use_dims = 2:3,
                     show_plots = F) +  theme_minimal(base_size = 8)
plotlist = list(current_points_plot, current_svd_plot)
cowplot::plot_grid(plotlist = plotlist, rel_widths = c(0.66, 0.33))

dso$plot_svd_history()

#Solution
dso$init_solution("random")
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  with_solution = TRUE,
  use_dims = list(2:3)
)

dso$default_optimization()

dso$plot_error_history()
dso$plot_negative_basis_change()
dso$plot_negative_proportions_change()

colnames(dso$st$proj$X) <- c("R1", "R2", "R3")
colnames(dso$st$proj$Omega) <- c("S1", "S2", "S3")
proj_solution_history <-
  dso$plot_projected(use_dims = (2:3),
                     wrap = F,
                     with_legend = F)
proj_solution_history[[1]] # genes in a space of samples
proj_solution_history[[2]] # samples in a space of genes

solution <- dso$finalize_solution()
result_W <- solution$W
result_H <- solution$H
true_H = sim_initial$proportions
ptp <- coerce_pred_true_props(solution$H, true_H)
plot_ptp_lines(ptp)

#M2
colors_v <- c("#1f77b4", "#ff7f0e", "#2ca02c")
colors <-
  which_marker(rownames(fData(dso$st$data)), dso$st$marker_genes)
plot_markers <-
  plot_projection_points(
    dso$st$proj,
    use_dims = (2:3),
    spaces = c("X"),
    pt_size = 1,
    color = colors
  ) +
  scale_color_manual(values = colors_v[1:3],
                     na.value = adjustcolor("grey70", alpha.f = 0.7)) +
  labs(col = "Marker Cell Type", x = "R2", y = "R3") +
  theme_bw(base_family = "sans", base_size = 12) +
  theme(
    legend.position = "right",
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

plot_markers