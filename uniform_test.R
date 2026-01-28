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

# Функция для добавления равномерного шума (предоставленная)
add_uniform_noise <- function(simulation, noise_range = c(-0.5, 0.5),
                              additive = TRUE,
                              protect_genes = c(), protect_samples = c()) {
  
  data <- simulation$data
  
  if (additive) {
    # Аддитивный равномерный шум
    noise <- matrix(
      stats::runif(n = length(data), min = noise_range[1], max = noise_range[2]),
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
    # Мультипликативный равномерный шум
    noise <- matrix(
      stats::runif(n = length(data), min = 1 + noise_range[1], max = 1 + noise_range[2]),
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
  simulation$uniform_noise_params <- list(
    range = noise_range,
    type = ifelse(additive, "additive", "multiplicative")
  )
  simulation$protected_genes <- protect_genes
  simulation$protected_samples <- protect_samples
  
  simulation
}

# Функция для оценки качества с равномерным шумом
evaluate_noise_effect <- function(range_vals, additive_vals, n_runs = 1) {
  results <- list()
  
  for (range_val in range_vals) {
    for (additive in additive_vals) {
      cat(sprintf("\n=== Testing: range = %.1f, additive = %s ===\n", 
                  range_val, ifelse(additive, "TRUE", "FALSE")))
      
      # Повторяем несколько раз для статистики
      run_metrics <- lapply(1:n_runs, function(run) {
        cat(sprintf("Run %d/%d... ", run, n_runs))
        
        # Создаем симуляцию
        set.seed(run)
        sim_initial <- create_simulation(n_genes = 10000, n_samples = 100, n_cell_types = 3)
        
        # Добавляем равномерный шум
        sim_noisy <- add_uniform_noise(
          sim_initial, 
          noise_range = c(-range_val, range_val), 
          additive = additive,
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
            range = range_val,
            additive = additive,
            run = run,
            pearson_mean = mean(metrics["pearson", ], na.rm = TRUE),
            rmse_mean = mean(metrics["rmse", ], na.rm = TRUE)
          )
          
          cat("OK\n")
          return(result_df)
          
        }, error = function(e) {
          cat(sprintf("ERROR: %s\n", e$message))
          return(data.frame(
            range = range_val,
            additive = additive,
            run = run,
            pearson_mean = NA,
            rmse_mean = NA
          ))
        })
      })
      
      results[[paste0("range", range_val, "_add", as.numeric(additive))]] <- do.call(rbind, run_metrics)
    }
  }
  
  # Объединяем все результаты
  final_results <- do.call(rbind, results)
  rownames(final_results) <- NULL
  return(final_results)
}

# Запускаем оценку с разными параметрами шума
range_vals <- c(0.1, 1, 5, 10, 20)
additive_vals <- c(TRUE, FALSE)

results <- evaluate_noise_effect(range_vals, additive_vals, n_runs = 1)

# Убираем NA значения для визуализации
results_clean <- na.omit(results)

# Визуализация результатов
library(ggplot2)

if (nrow(results_clean) > 0) {
  # Преобразуем additive в фактор для лучшей визуализации
  results_clean$additive <- factor(results_clean$additive, 
                                   levels = c(TRUE, FALSE),
                                   labels = c("Additive", "Multiplicative"))
  
  # График для Pearson correlation
  p1 <- ggplot(results_clean, aes(x = factor(range), y = pearson_mean, 
                                  color = additive, group = additive)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Range (±value)", y = "Pearson Correlation", 
         color = "Noise Type", 
         title = "Влияние равномерного шума на точность деконволюции (Pearson)") +
    theme_minimal() +
    ylim(-0.1, 1.1)
  
  # График для RMSE
  p2 <- ggplot(results_clean, aes(x = factor(range), y = rmse_mean,
                                  color = additive, group = additive)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Range (±value)", y = "RMSE", 
         color = "Noise Type",
         title = "Влияние равномерного шума на точность деконволюции (RMSE)") +
    theme_minimal()
  
  # Показываем графики
  print(p1)
  print(p2)
  
  # Heatmap средних значений
  library(reshape2)
  
  # Агрегируем данные
  agg_data <- aggregate(cbind(pearson_mean, rmse_mean) ~ range + additive, 
                        data = results_clean, FUN = mean, na.action = na.pass)
  
  # Heatmap для Pearson
  pearson_mat <- acast(agg_data, range ~ additive, 
                       value.var = "pearson_mean")
  
  if (!all(is.na(pearson_mat))) {
    heatmap_pearson <- ggplot(melt(pearson_mat), 
                              aes(x = factor(Var2), y = factor(Var1), fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                           midpoint = 0.5, limits = c(-0.1, 1), na.value = "grey90") +
      labs(x = "Noise Type", y = "Range (±value)", 
           fill = "Pearson",
           title = "Средняя корреляция Пирсона") +
      theme_minimal()
    
    print(heatmap_pearson)
  }
  
  # Heatmap для RMSE
  rmse_mat <- acast(agg_data, range ~ additive, 
                    value.var = "rmse_mean")
  
  if (!all(is.na(rmse_mat))) {
    heatmap_rmse <- ggplot(melt(rmse_mat), 
                           aes(x = factor(Var2), y = factor(Var1), fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 3)), color = "black", size = 4) +
      scale_fill_gradient(low = "green", high = "red", na.value = "grey90") +
      labs(x = "Noise Type", y = "Range (±value)", 
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