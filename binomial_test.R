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

# Модифицированная функция для добавления шума с защитой от нулевых строк
add_negative_binomial_noise_safe <- function(simulation, dispersion = 0.1,
                                             scaling_factor = 1,
                                             min_nonzero = 1,  # минимальное ненулевое значение в строке
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
  # Защита от нулевых строк: если вся строка нулевая, добавляем небольшой шум
  zero_rows <- which(rowSums(noisy_data) == 0)
  if (length(zero_rows) > 0) {
    cat(sprintf("Found %d zero rows, adding small noise...\n", length(zero_rows)))
    for (row in zero_rows) {
      noisy_data[row, ] <- stats::rpois(ncol(noisy_data), lambda = min_nonzero)
    }
  }
  
  simulation$data <- noisy_data
  simulation$nb_noise_params <- list(
    dispersion = dispersion,
    scaling_factor = scaling_factor,
    size = size_param,
    fixed_zero_rows = length(zero_rows)
  )
  
  simulation
}

# Функция для оценки качества
evaluate_noise_effect <- function(dispersion_vals, scaling_vals, n_runs = 3) {
  results <- list()
  
  for (disp in dispersion_vals) {
    for (scale in scaling_vals) {
      cat(sprintf("\n=== Testing: dispersion = %.1f, scaling = %.1f ===\n", disp, scale))
      
      # Повторяем несколько раз для статистики
      run_metrics <- lapply(1:n_runs, function(run) {
        cat(sprintf("Run %d/%d... ", run, n_runs))
        
        # Создаем симуляцию с меньшим числом генов для скорости
        set.seed(run)
        sim_initial <- create_simulation(n_genes = 10000, n_samples = 100, n_cell_types = 3)
        
        # Добавляем шум с защитой
        sim_noisy <- add_negative_binomial_noise_safe(
          sim_initial, 
          dispersion = disp, 
          scaling_factor = scale, 
          min_nonzero = 0.1,
          seed = run
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
            dispersion = disp,
            scaling_factor = scale,
            run = run,
            pearson_mean = mean(metrics["pearson", ], na.rm = TRUE),
            rmse_mean = mean(metrics["rmse", ], na.rm = TRUE)
          )
          
          cat("OK\n")
          return(result_df)
          
        }, error = function(e) {
          cat(sprintf("ERROR: %s\n", e$message))
          return(data.frame(
            dispersion = disp,
            scaling_factor = scale,
            run = run,
            pearson_mean = NA,
            rmse_mean = NA
          ))
        })
      })
      
      results[[paste0("disp", disp, "_scale", scale)]] <- do.call(rbind, run_metrics)
    }
  }
  
  # Объединяем все результаты
  final_results <- do.call(rbind, results)
  rownames(final_results) <- NULL
  return(final_results)
}

# Запускаем оценку с разными параметрами шума
dispersion_vals <- c(0.1, 0.5, 1, 10, 100)
scaling_vals <- c(0.1, 0.5, 1, 5, 10)

results <- evaluate_noise_effect(dispersion_vals, scaling_vals, n_runs = 1)

# Убираем NA значения для визуализации
results_clean <- na.omit(results)

# Визуализация результатов
library(ggplot2)

if (nrow(results_clean) > 0) {
  # График для Pearson correlation
  p1 <- ggplot(results_clean, aes(x = factor(dispersion), y = pearson_mean, 
                                  color = factor(scaling_factor), group = scaling_factor)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Dispersion (φ)", y = "Pearson Correlation", 
         color = "Scaling Factor", 
         title = "Влияние шума на точность деконволюции (Pearson)") +
    theme_minimal() +
    ylim(-0.1, 1.1)
  
  # График для RMSE
  p2 <- ggplot(results_clean, aes(x = factor(dispersion), y = rmse_mean,
                                  color = factor(scaling_factor), group = scaling_factor)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Dispersion (φ)", y = "RMSE", 
         color = "Scaling Factor",
         title = "Влияние шума на точность деконволюции (RMSE)") +
    theme_minimal()
  
  # Показываем графики
  print(p1)
  print(p2)
  
  # Heatmap средних значений
  library(reshape2)
  
  # Агрегируем данные
  agg_data <- aggregate(cbind(pearson_mean, rmse_mean) ~ dispersion + scaling_factor, 
                        data = results_clean, FUN = mean, na.action = na.pass)
  
  # Heatmap для Pearson
  pearson_mat <- acast(agg_data, dispersion ~ scaling_factor, 
                       value.var = "pearson_mean")
  
  if (!all(is.na(pearson_mat))) {
    heatmap_pearson <- ggplot(melt(pearson_mat), 
                              aes(x = factor(Var2), y = factor(Var1), fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                           midpoint = 0.5, limits = c(-0.1, 1), na.value = "grey90") +
      labs(x = "Scaling Factor", y = "Dispersion (φ)", 
           fill = "Pearson",
           title = "Средняя корреляция Пирсона") +
      theme_minimal()
    
    print(heatmap_pearson)
  }
  
  # Heatmap для RMSE
  rmse_mat <- acast(agg_data, dispersion ~ scaling_factor, 
                    value.var = "rmse_mean")
  
  if (!all(is.na(rmse_mat))) {
    heatmap_rmse <- ggplot(melt(rmse_mat), 
                           aes(x = factor(Var2), y = factor(Var1), fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 3)), color = "black", size = 4) +
      scale_fill_gradient(low = "green", high = "red", na.value = "grey90") +
      labs(x = "Scaling Factor", y = "Dispersion (φ)", 
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
