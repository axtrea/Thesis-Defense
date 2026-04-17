library(DualSimplex)
library(dplyr)
library(ggplot2)
library(matrixStats)

# ===============================
# Noise function
# ===============================
add_gaussian_noise <- function(simulation, noise_deviation, additive = FALSE) {
  data <- simulation$data
  
  noise <- matrix(
    stats::rnorm(n = length(data), mean = 0, sd = noise_deviation),
    nrow = nrow(data),
    ncol = ncol(data)
  )
  noise <- 2^noise
  
  if (additive) {
    noisy_data <- data + noise
  } else {
    noisy_data <- data * noise
  }
  
  noisy_data[noisy_data < 0] <- 0
  simulation$data <- noisy_data
  simulation
}

# ===============================
# Metrics
# ===============================
compute_metrics <- function(pred, true) {
  rmse <- sqrt(mean((pred - true)^2))
  corrs <- sapply(1:nrow(true), function(i)
    cor(pred[i, ], true[i, ]))
  pearson <- mean(corrs)
  
  list(RMSE = rmse, Pearson = pearson)
}

# ===============================
# Experiment settings
# ===============================
noise_levels <- c(3, 4, 4.5, 5, 6)
n_repeats <- 2

results <- data.frame()

# ===============================
# Main experiment loop
# ===============================
for (noise in noise_levels) {
  for (rep in 1:n_repeats) {
    
    cat("Noise:", noise, " Rep:", rep, "\n")
    
    # --- Simulation ---
    sim <- create_simulation(
      n_genes = 10000,
      n_samples = 100,
      n_cell_types = 3,
      with_marker_genes = FALSE
    )
    
    sim <- add_gaussian_noise(sim, noise_deviation = noise, additive = TRUE)
    
    data_raw <- sim$data
    true_props <- sim$proportions
    
    # ===============================
    # DualSimplex с фильтром distance_filter
    # ===============================
    dso_dist <- DualSimplexSolver$new()
    dso_dist$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
    dso_dist$project(3)
    
    # Применяем distance_filter (порог 0.10)
    dso_dist$distance_filter(plane_d_lt = 0.15, genes = TRUE)
    
    dso_dist$project(3)
    dso_dist$init_solution("random")
    dso_dist$default_optimization()
    solution_dist <- dso_dist$finalize_solution()
    
    # Проверка и сохранение результатов для distance_filter
    if (!is.null(solution_dist$H) && !any(is.na(solution_dist$H))) {
      ptp_ds <- coerce_pred_true_props(solution_dist$H, true_props)
      metrics_ds <- compute_metrics(ptp_ds[[2]], ptp_ds[[1]])
      
      results <- rbind(results, data.frame(
        Method = "DS_distance",
        Noise = noise,
        RMSE = metrics_ds$RMSE,
        Pearson = metrics_ds$Pearson
      ))
    } else {
      cat("DS_distance: пропуск повтора (NA или NULL) при noise =", noise, "rep =", rep, "\n")
    }
    
    # ===============================
    # DualSimplex с фильтром basic_filter (MAD)
    # ===============================
    dso_basic <- DualSimplexSolver$new()
    dso_basic$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
    dso_basic$project(3)
    
    # Применяем basic_filter с порогом log_mad_gt = 3
    tryCatch({
      dso_basic$basic_filter(log_mad_gt = 3, genes = TRUE)
    }, error = function(e) {
      cat("Ошибка в basic_filter: ", e$message, "\n")
      # В случае ошибки создадим пустое решение, чтобы пропустить повтор
      solution_basic <<- NULL
    })
    
    # Если фильтрация прошла успешно, продолжаем
    if (exists("solution_basic") && is.null(solution_basic)) {
      # ошибка уже обработана, пропускаем
      cat("DS_basic: пропуск повтора из-за ошибки фильтрации при noise =", noise, "rep =", rep, "\n")
      next
    }
    
    dso_basic$project(3)
    dso_basic$init_solution("random")
    dso_basic$default_optimization()
    solution_basic <- dso_basic$finalize_solution()
    
    if (!is.null(solution_basic$H) && !any(is.na(solution_basic$H))) {
      ptp_ds <- coerce_pred_true_props(solution_basic$H, true_props)
      metrics_ds <- compute_metrics(ptp_ds[[2]], ptp_ds[[1]])
      
      results <- rbind(results, data.frame(
        Method = "DS_basic",
        Noise = noise,
        RMSE = metrics_ds$RMSE,
        Pearson = metrics_ds$Pearson
      ))
    } else {
      cat("DS_basic: пропуск повтора (NA или NULL) при noise =", noise, "rep =", rep, "\n")
    }
    
  } # конец цикла по повторениям (rep)
} # конец цикла по уровням шума (noise)

# ===============================
# Фильтрация выбросов (RMSE > 2) - опционально
# ===============================
cat("До фильтрации:", nrow(results), "строк\n")
results <- results %>% filter(RMSE <= 2)
cat("После фильтрации:", nrow(results), "строк\n")

# ===============================
# Aggregate results
# ===============================
summary_results <- results %>%
  group_by(Method, Noise) %>%
  summarise(
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    RMSE_sd = sd(RMSE, na.rm = TRUE),
    Pearson_mean = mean(Pearson, na.rm = TRUE),
    Pearson_sd = sd(Pearson, na.rm = TRUE),
    .groups = "drop"
  )

# ===============================
# Построение графиков
# ===============================

library(ggplot2)

# Цвета для двух методов
method_colors <- c("DS_distance" = "steelblue", "DS_basic" = "coral")

# График для RMSE
p_rmse <- ggplot(summary_results, aes(x = Noise, y = RMSE_mean, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = RMSE_mean - RMSE_sd, ymax = RMSE_mean + RMSE_sd), width = 0.2) +
  labs(
    title = "Сравнение фильтров DualSimplex по RMSE",
    x = "Уровень шума (noise)",
    y = "RMSE (среднее ± SD)"
  ) +
  theme_minimal() +
  scale_color_manual(values = method_colors) +
  theme(legend.position = "bottom")

# График для Pearson correlation
p_pearson <- ggplot(summary_results, aes(x = Noise, y = Pearson_mean, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Pearson_mean - Pearson_sd, ymax = Pearson_mean + Pearson_sd), width = 0.2) +
  labs(
    title = "Сравнение фильтров DualSimplex по коэффициенту корреляции Пирсона",
    x = "Уровень шума (noise)",
    y = "Pearson correlation (среднее ± SD)"
  ) +
  theme_minimal() +
  scale_color_manual(values = method_colors) +
  theme(legend.position = "bottom")

# Вывод графиков
print(p_rmse)
print(p_pearson)

# ===============================
# Boxplot для RMSE
# ===============================
p_rmse_box <- ggplot(results, aes(x = factor(Noise), y = RMSE, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             alpha = 0.5, size = 1) +
  labs(
    title = "Сравнение фильтров DualSimplex по RMSE",
    x = "Уровень шума",
    y = "RMSE"
  ) +
  theme_minimal() +
  scale_fill_manual(values = method_colors) +
  theme(legend.position = "bottom")

# ===============================
# Boxplot для Pearson correlation
# ===============================
p_pearson_box <- ggplot(results, aes(x = factor(Noise), y = Pearson, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             alpha = 0.5, size = 1) +
  labs(
    title = "Сравнение фильтров DualSimplex по коэффициенту корреляции Пирсона",
    x = "Уровень шума",
    y = "Pearson correlation"
  ) +
  theme_minimal() +
  scale_fill_manual(values = method_colors) +
  theme(legend.position = "bottom")

# Вывод графиков
print(p_rmse_box)
print(p_pearson_box)