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

# Функция для добавления шума "соль и перец" (предоставленная)
add_salt_pepper_noise <- function(simulation, probability = 0.01, 
                                  salt_value = NULL, pepper_value = NULL,
                                  salt_ratio = 0.5,
                                  protect_genes = c(), protect_samples = c()) {
  
  data <- simulation$data
  
  # Определяем значения по умолчанию
  if (is.null(salt_value)) {
    salt_value <- quantile(data, 0.95) * 2  # Высокое значение
  }
  
  if (is.null(pepper_value)) {
    pepper_value <- quantile(data, 0.05) / 2  # Низкое значение
  }
  
  # Создаем маску шума
  noise_mask <- matrix(
    stats::runif(n = length(data)),
    nrow = nrow(data),
    ncol = ncol(data)
  )
  
  # Применяем шум
  noisy_data <- data
  
  # Находим позиции для шума
  salt_positions <- noise_mask < (probability * salt_ratio)
  pepper_positions <- (noise_mask >= (probability * salt_ratio)) & 
    (noise_mask < probability)
  
  # Защита генов и образцов
  if (length(protect_genes) > 0) {
    salt_positions[protect_genes, ] <- FALSE
    pepper_positions[protect_genes, ] <- FALSE
  }
  
  if (length(protect_samples) > 0) {
    salt_positions[, protect_samples] <- FALSE
    pepper_positions[, protect_samples] <- FALSE
  }
  
  # Применяем "соль" и "перец"
  noisy_data[salt_positions] <- salt_value
  noisy_data[pepper_positions] <- pepper_value
  
  simulation$data <- noisy_data
  simulation$salt_pepper_noise_params <- list(
    probability = probability,
    salt_value = salt_value,
    pepper_value = pepper_value,
    salt_ratio = salt_ratio
  )
  simulation$protected_genes <- protect_genes
  simulation$protected_samples <- protect_samples
  
  simulation
}

# Функция для оценки качества с шумом "соль и перец"
evaluate_noise_effect <- function(probability_vals, salt_ratio_vals, n_runs = 3) {
  results <- list()
  
  for (prob in probability_vals) {
    for (ratio in salt_ratio_vals) {
      cat(sprintf("\n=== Testing: probability = %.3f, salt_ratio = %.1f ===\n", 
                  prob, ratio))
      
      # Повторяем несколько раз для статистики
      run_metrics <- lapply(1:n_runs, function(run) {
        cat(sprintf("Run %d/%d... ", run, n_runs))
        
        # Создаем симуляцию с меньшим числом генов для скорости
        set.seed(run)
        sim_initial <- create_simulation(n_genes = 10000, n_samples = 100, n_cell_types = 3)
        
        # Добавляем шум "соль и перец"
        sim_noisy <- add_salt_pepper_noise(
          sim_initial, 
          probability = prob, 
          salt_ratio = ratio,
          protect_genes = c(),
          protect_samples = c()
        )
        
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
            probability = prob,
            salt_ratio = ratio,
            run = run,
            pearson_mean = mean(metrics["pearson", ], na.rm = TRUE),
            rmse_mean = mean(metrics["rmse", ], na.rm = TRUE)
          )
          
          cat("OK\n")
          return(result_df)
          
        }, error = function(e) {
          cat(sprintf("ERROR: %s\n", e$message))
          return(data.frame(
            probability = prob,
            salt_ratio = ratio,
            run = run,
            pearson_mean = NA,
            rmse_mean = NA
          ))
        })
      })
      
      results[[paste0("prob", prob, "_ratio", ratio)]] <- do.call(rbind, run_metrics)
    }
  }
  
  # Объединяем все результаты
  final_results <- do.call(rbind, results)
  rownames(final_results) <- NULL
  return(final_results)
}

# Запускаем оценку с разными параметрами шума
probability_vals <- c(0.01, 0.05, 0.1, 0.2, 0.3)
salt_ratio_vals <- c(0.2, 0.5, 0.8)

results <- evaluate_noise_effect(probability_vals, salt_ratio_vals, n_runs = 2)

# Убираем NA значения для визуализации
results_clean <- na.omit(results)

# Визуализация результатов
library(ggplot2)

if (nrow(results_clean) > 0) {
  # График для Pearson correlation
  p1 <- ggplot(results_clean, aes(x = factor(probability), y = pearson_mean, 
                                  color = factor(salt_ratio), group = salt_ratio)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Probability", y = "Pearson Correlation", 
         color = "Salt Ratio", 
         title = "Влияние шума 'соль и перец' на точность деконволюции (Pearson)") +
    theme_minimal() +
    ylim(-0.1, 1.1)
  
  # График для RMSE
  p2 <- ggplot(results_clean, aes(x = factor(probability), y = rmse_mean,
                                  color = factor(salt_ratio), group = salt_ratio)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Probability", y = "RMSE", 
         color = "Salt Ratio",
         title = "Влияние шума 'соль и перец' на точность деконволюции (RMSE)") +
    theme_minimal()
  
  # Показываем графики
  print(p1)
  print(p2)
  
  # Heatmap средних значений
  library(reshape2)
  
  # Агрегируем данные
  agg_data <- aggregate(cbind(pearson_mean, rmse_mean) ~ probability + salt_ratio, 
                        data = results_clean, FUN = mean, na.action = na.pass)
  
  # Heatmap для Pearson
  pearson_mat <- acast(agg_data, probability ~ salt_ratio, 
                       value.var = "pearson_mean")
  
  if (!all(is.na(pearson_mat))) {
    heatmap_pearson <- ggplot(melt(pearson_mat), 
                              aes(x = factor(Var2), y = factor(Var1), fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                           midpoint = 0.5, limits = c(-0.1, 1), na.value = "grey90") +
      labs(x = "Salt Ratio", y = "Probability", 
           fill = "Pearson",
           title = "Средняя корреляция Пирсона") +
      theme_minimal()
    
    print(heatmap_pearson)
  }
  
  # Heatmap для RMSE
  rmse_mat <- acast(agg_data, probability ~ salt_ratio, 
                    value.var = "rmse_mean")
  
  if (!all(is.na(rmse_mat))) {
    heatmap_rmse <- ggplot(melt(rmse_mat), 
                           aes(x = factor(Var2), y = factor(Var1), fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 3)), color = "black", size = 4) +
      scale_fill_gradient(low = "green", high = "red", na.value = "grey90") +
      labs(x = "Salt Ratio", y = "Probability", 
           fill = "RMSE",
           title = "Средняя RMSE") +
      theme_minimal()
    
    print(heatmap_rmse)
  }
  
  # Сводная статистика
  cat("\n=== Сводная статистика ===\n")
  print(agg_data)
  
} else {
  cat("\nВсе запуски завершились с ошибками. Попробуйте уменьшить параметры шума.\n")
}