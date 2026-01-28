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

add_poisson_noise <- function(simulation, scaling_factor = 1,
                              protect_genes = c(), protect_samples = c()) {
  
  data <- simulation$data
  
  # Для Пуассоновского шума добавляем шум, пропорциональный корню из значения
  # Это аппроксимация для больших значений
  poisson_noise <- matrix(
    stats::rnorm(n = length(data), mean = 0, sd = 1),
    nrow = nrow(data),
    ncol = ncol(data)
  )
  
  # Масштабируем шум в соответствии с корнем из значения
  scaled_noise <- poisson_noise * sqrt(abs(data)) * scaling_factor
  
  
  # Защита генов и образцов
  if (length(protect_genes) > 0) {
    scaled_noise[protect_genes, ] <- 0
  }
  if (length(protect_samples) > 0) {
    scaled_noise[, protect_samples] <- 0
  }
  
  noisy_data <- data + scaled_noise
  
  # Обеспечение неотрицательности
  noisy_data[noisy_data < 0] <- 0
  
  simulation$data <- noisy_data
  simulation$poisson_noise_params <- list(scaling_factor = scaling_factor)
  simulation$protected_genes <- protect_genes
  simulation$protected_samples <- protect_samples
  
  simulation
}

# Функция для оценки качества с Пуассоновским шумом
evaluate_poisson_noise_effect <- function(scaling_vals, n_runs = 3) {
  results <- list()
  
  for (scale in scaling_vals) {
    cat(sprintf("\n=== Testing: scaling factor = %.1f ===\n", scale))
    
    # Повторяем несколько раз для статистики
    run_metrics <- lapply(1:n_runs, function(run) {
      cat(sprintf("Run %d/%d... ", run, n_runs))
      
      # Создаем симуляцию с меньшим числом генов для скорости
      set.seed(run)
      sim_initial <- create_simulation(n_genes = 10000, n_samples = 100, n_cell_types = 3)
      
      # Добавляем Пуассоновский шум
      sim_noisy <- add_poisson_noise(
        sim_initial, 
        scaling_factor = scale, 
        protect_genes = c(),
        protect_samples = c()
      )
      
      # Проверяем нулевые строки
      zero_rows_final <- sum(rowSums(sim_noisy$data) == 0)
      if (zero_rows_final > 0) {
        warning(sprintf("Still %d zero rows after noise addition!", zero_rows_final))
      }
      
      tryCatch({
        # Запускаем DualSimplex
        dso <- DualSimplexSolver$new()
        dso$set_data(sim_noisy$data, max_dim = 10)
        dso$project(3)
        dso$init_solution("random")
        dso$default_optimization()
        solution <- dso$finalize_solution()
        
        # Сравниваем с истинными пропорциями
        ptp <- coerce_pred_true_props(solution$H, sim_initial$proportions)
        pred_H <- ptp[[1]]
        true_H <- ptp[[2]]
        
        # Вычисляем метрики для каждого типа клеток
        metrics <- sapply(1:nrow(pred_H), function(i) {
          pred <- pred_H[i, ]
          true <- true_H[i, ]
          c(
            pearson = ifelse(sd(pred) > 0 & sd(true) > 0, 
                             cor(pred, true, method = "pearson"), 
                             NA),
            rmse = sqrt(mean((pred - true)^2))
          )
        })
        
        # Усредняем по типам клеток
        result_df <- data.frame(
          scaling_factor = scale,
          run = run,
          pearson_mean = mean(metrics["pearson", ], na.rm = TRUE),
          rmse_mean = mean(metrics["rmse", ], na.rm = TRUE),
          zero_rows = zero_rows_final
        )
        
        cat("OK\n")
        return(result_df)
        
      }, error = function(e) {
        cat(sprintf("ERROR: %s\n", e$message))
        return(data.frame(
          scaling_factor = scale,
          run = run,
          pearson_mean = NA,
          rmse_mean = NA,
          zero_rows = zero_rows_final
        ))
      })
    })
    
    results[[paste0("scale", scale)]] <- do.call(rbind, run_metrics)
  }
  
  # Объединяем все результаты
  final_results <- do.call(rbind, results)
  rownames(final_results) <- NULL
  return(final_results)
}

# Запускаем оценку с разными параметрами Пуассоновского шума
scaling_vals <- c(0.1, 0.5, 1, 5, 10, 100)  # Коэффициенты масштабирования для Пуассона

results <- evaluate_poisson_noise_effect(scaling_vals, n_runs = 2)

# Убираем NA значения для визуализации
results_clean <- na.omit(results)

# Визуализация результатов
if (nrow(results_clean) > 0) {
  # График для Pearson correlation
  p1 <- ggplot(results_clean, aes(x = factor(scaling_factor), y = pearson_mean)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red", size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
    labs(x = "Scaling Factor (λ multiplier)", y = "Pearson Correlation", 
         title = "Влияние Пуассоновского шума на точность деконволюции (Pearson)",
         subtitle = "Для Пуассона: Var = μ, Scaling Factor умножает μ") +
    theme_minimal() +
    ylim(-0.1, 1.1)
  
  # График для RMSE
  p2 <- ggplot(results_clean, aes(x = factor(scaling_factor), y = rmse_mean)) +
    geom_boxplot(fill = "lightcoral", alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "blue", size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "blue") +
    labs(x = "Scaling Factor (λ multiplier)", y = "RMSE",
         title = "Влияние Пуассоновского шума на точность деконволюции (RMSE)") +
    theme_minimal()
  
  # График для количества нулевых строк
  p3 <- ggplot(results_clean, aes(x = factor(scaling_factor), y = zero_rows)) +
    geom_boxplot(fill = "lightgreen", alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "darkgreen", size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "darkgreen") +
    labs(x = "Scaling Factor (λ multiplier)", y = "Количество нулевых строк",
         title = "Влияние масштабирования на количество нулевых строк") +
    theme_minimal()
  
  # Показываем графики
  print(p1)
  print(p2)
  print(p3)
  
  # Heatmap средних значений
  library(reshape2)
  
  # Агрегируем данные
  agg_data <- aggregate(cbind(pearson_mean, rmse_mean, zero_rows) ~ scaling_factor, 
                        data = results_clean, FUN = mean, na.action = na.pass)
  
  # Heatmap для Pearson
  pearson_mat <- matrix(agg_data$pearson_mean, nrow = 1)
  colnames(pearson_mat) <- agg_data$scaling_factor
  
  if (!all(is.na(pearson_mat))) {
    pearson_melt <- melt(pearson_mat)
    
    heatmap_pearson <- ggplot(pearson_melt, 
                              aes(x = factor(Var2), y = 1, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 3)), color = "black", size = 4) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                           midpoint = 0.5, limits = c(-0.1, 1), na.value = "grey90") +
      labs(x = "Scaling Factor", y = "", 
           fill = "Pearson",
           title = "Средняя корреляция Пирсона для Пуассоновского шума") +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    
    print(heatmap_pearson)
  }
  
  # Heatmap для RMSE
  rmse_mat <- matrix(agg_data$rmse_mean, nrow = 1)
  colnames(rmse_mat) <- agg_data$scaling_factor
  
  if (!all(is.na(rmse_mat))) {
    rmse_melt <- melt(rmse_mat)
    
    heatmap_rmse <- ggplot(rmse_melt, 
                           aes(x = factor(Var2), y = 1, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 3)), color = "black", size = 4) +
      scale_fill_gradient(low = "green", high = "red", na.value = "grey90") +
      labs(x = "Scaling Factor", y = "", 
           fill = "RMSE",
           title = "Средняя RMSE для Пуассоновского шума") +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    
    print(heatmap_rmse)
  }
  
  # Сводная статистика
  cat("\n=== Сводная статистика ===\n")
  print(agg_data)
  
  # Дополнительный анализ: тренды
  cat("\n=== Анализ трендов ===\n")
  
  # Линейная регрессия для анализа трендов
  if (length(unique(results_clean$scaling_factor)) > 2) {
    # Для Pearson
    lm_pearson <- lm(pearson_mean ~ scaling_factor, data = results_clean)
    cat("Линейная регрессия для Pearson:\n")
    print(summary(lm_pearson))
    
    # Для RMSE
    lm_rmse <- lm(rmse_mean ~ scaling_factor, data = results_clean)
    cat("\nЛинейная регрессия для RMSE:\n")
    print(summary(lm_rmse))
    
    # График с линией регрессии
    p_reg_pearson <- ggplot(results_clean, aes(x = scaling_factor, y = pearson_mean)) +
      geom_point(alpha = 0.5, size = 3) +
      geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
      labs(x = "Scaling Factor", y = "Pearson Correlation",
           title = "Линейная регрессия: Pearson vs Scaling Factor") +
      theme_minimal()
    
    p_reg_rmse <- ggplot(results_clean, aes(x = scaling_factor, y = rmse_mean)) +
      geom_point(alpha = 0.5, size = 3) +
      geom_smooth(method = "lm", se = TRUE, color = "red", fill = "lightcoral") +
      labs(x = "Scaling Factor", y = "RMSE",
           title = "Линейная регрессия: RMSE vs Scaling Factor") +
      theme_minimal()
    
    print(p_reg_pearson)
    print(p_reg_rmse)
  }
  
} else {
  cat("\nВсе запуски завершились с ошибками. Попробуйте изменить параметры шума.\n")
}
