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

# Функция для добавления dropout шума (предоставленная)
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

# Функция для оценки качества с dropout шумом
evaluate_noise_effect <- function(dropout_rates, mean_based_vals, n_runs = 3) {
  results <- list()
  
  for (rate in dropout_rates) {
    for (mean_based in mean_based_vals) {
      cat(sprintf("\n=== Testing: dropout_rate = %.1f, mean_based = %s ===\n", 
                  rate, ifelse(mean_based, "TRUE", "FALSE")))
      
      # Повторяем несколько раз для статистики
      run_metrics <- lapply(1:n_runs, function(run) {
        cat(sprintf("Run %d/%d... ", run, n_runs))
        
        # Создаем симуляцию с меньшим числом генов для скорости
        set.seed(run)
        sim_initial <- create_simulation(n_genes = 10000, n_samples = 100, n_cell_types = 3)
        
        # Добавляем dropout шум
        sim_noisy <- add_dropout_noise(
          sim_initial, 
          dropout_rate = rate,
          mean_based = mean_based,
          protect_genes = c(),
          protect_samples = c(),
          seed = run
        )
        
        # Проверяем количество нулей
        zeros_proportion <- sum(sim_noisy$data == 0) / length(sim_noisy$data)
        cat(sprintf("Zeros: %.1f%% ", zeros_proportion * 100))
        
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
            dropout_rate = rate,
            mean_based = mean_based,
            run = run,
            zeros_proportion = zeros_proportion,
            pearson_mean = mean(metrics["pearson", ], na.rm = TRUE),
            rmse_mean = mean(metrics["rmse", ], na.rm = TRUE)
          )
          
          cat("OK\n")
          return(result_df)
          
        }, error = function(e) {
          cat(sprintf("ERROR: %s\n", e$message))
          return(data.frame(
            dropout_rate = rate,
            mean_based = mean_based,
            run = run,
            zeros_proportion = zeros_proportion,
            pearson_mean = NA,
            rmse_mean = NA
          ))
        })
      })
      
      results[[paste0("rate", rate, "_mean", as.numeric(mean_based))]] <- do.call(rbind, run_metrics)
    }
  }
  
  # Объединяем все результаты
  final_results <- do.call(rbind, results)
  rownames(final_results) <- NULL
  return(final_results)
}

# Запускаем оценку с разными параметрами шума
dropout_rates <- c(0.1, 0.3, 0.5, 0.7, 0.9)
mean_based_vals <- c(TRUE, FALSE)

results <- evaluate_noise_effect(dropout_rates, mean_based_vals, n_runs = 2)

# Убираем NA значения для визуализации
results_clean <- na.omit(results)

# Визуализация результатов
library(ggplot2)

if (nrow(results_clean) > 0) {
  # Преобразуем mean_based в фактор для лучшей визуализации
  results_clean$mean_based <- factor(results_clean$mean_based, 
                                     levels = c(TRUE, FALSE),
                                     labels = c("Mean-based", "Uniform"))
  
  # График для Pearson correlation
  p1 <- ggplot(results_clean, aes(x = factor(dropout_rate), y = pearson_mean, 
                                  color = mean_based, group = mean_based)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Dropout Rate", y = "Pearson Correlation", 
         color = "Dropout Type", 
         title = "Влияние dropout шума на точность деконволюции (Pearson)") +
    theme_minimal() +
    ylim(-0.1, 1.1)
  
  # График для RMSE
  p2 <- ggplot(results_clean, aes(x = factor(dropout_rate), y = rmse_mean,
                                  color = mean_based, group = mean_based)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Dropout Rate", y = "RMSE", 
         color = "Dropout Type",
         title = "Влияние dropout шума на точность деконволюции (RMSE)") +
    theme_minimal()
  
  # График для количества нулей
  p3 <- ggplot(results_clean, aes(x = factor(dropout_rate), y = zeros_proportion,
                                  color = mean_based, group = mean_based)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Dropout Rate", y = "Proportion of Zeros", 
         color = "Dropout Type",
         title = "Фактическая доля нулей в данных") +
    theme_minimal() +
    scale_y_continuous(labels = scales::percent)
  
  # Показываем графики
  print(p1)
  print(p2)
  print(p3)
  
  # Heatmap средних значений
  library(reshape2)
  
  # Агрегируем данные
  agg_data <- aggregate(cbind(pearson_mean, rmse_mean, zeros_proportion) ~ dropout_rate + mean_based, 
                        data = results_clean, FUN = mean, na.action = na.pass)
  
  # Heatmap для Pearson
  pearson_mat <- acast(agg_data, dropout_rate ~ mean_based, 
                       value.var = "pearson_mean")
  
  if (!all(is.na(pearson_mat))) {
    heatmap_pearson <- ggplot(melt(pearson_mat), 
                              aes(x = factor(Var2), y = factor(Var1), fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                           midpoint = 0.5, limits = c(-0.1, 1), na.value = "grey90") +
      labs(x = "Dropout Type", y = "Dropout Rate", 
           fill = "Pearson",
           title = "Средняя корреляция Пирсона") +
      theme_minimal()
    
    print(heatmap_pearson)
  }
  
  # Heatmap для RMSE
  rmse_mat <- acast(agg_data, dropout_rate ~ mean_based, 
                    value.var = "rmse_mean")
  
  if (!all(is.na(rmse_mat))) {
    heatmap_rmse <- ggplot(melt(rmse_mat), 
                           aes(x = factor(Var2), y = factor(Var1), fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(value, 3)), color = "black", size = 4) +
      scale_fill_gradient(low = "green", high = "red", na.value = "grey90") +
      labs(x = "Dropout Type", y = "Dropout Rate", 
           fill = "RMSE",
           title = "Средняя RMSE") +
      theme_minimal()
    
    print(heatmap_rmse)
  }
  
  # Сводная статистика
  cat("\n=== Сводная статистика ===\n")
  print(agg_data)
  
  # Дополнительный анализ: соотношение между dropout_rate и фактическим количеством нулей
  cat("\n=== Анализ эффективности dropout модели ===\n")
  for (type in levels(results_clean$mean_based)) {
    type_data <- results_clean[results_clean$mean_based == type, ]
    lm_model <- lm(zeros_proportion ~ dropout_rate, data = type_data)
    cat(sprintf("\nДля %s dropout:\n", type))
    cat(sprintf("  R-squared: %.3f\n", summary(lm_model)$r.squared))
    cat(sprintf("  Фактические нули = %.3f + %.3f * dropout_rate\n", 
                coef(lm_model)[1], coef(lm_model)[2]))
  }
  
} else {
  cat("\nВсе запуски завершились с ошибками. Попробуйте уменьшить параметры шума.\n")
}

# Дополнительная функция для визуализации эффекта dropout
visualize_dropout_effect <- function() {
  # Создаем симуляцию
  sim_initial <- create_simulation(n_genes = 1000, n_samples = 50, n_cell_types = 3)
  
  # Применяем разные уровни dropout
  dropout_levels <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  
  par(mfrow = c(2, 3))
  
  # Исходные данные
  hist(sim_initial$data, main = "Исходные данные", xlab = "Expression", breaks = 30, 
       xlim = c(0, max(sim_initial$data)))
  
  for (rate in dropout_levels) {
    sim_dropout <- add_dropout_noise(sim_initial, dropout_rate = rate, mean_based = TRUE, seed = 123)
    
    # Вычисляем статистику
    original_nonzero <- sum(sim_initial$data > 0)
    dropout_nonzero <- sum(sim_dropout$data > 0)
    zeros_proportion <- sum(sim_dropout$data == 0) / length(sim_dropout$data)
    
    # Гистограмма
    hist(sim_dropout$data, 
         main = sprintf("Dropout rate = %.1f\nNon-zero: %d → %d (%.1f%% zeros)", 
                        rate, original_nonzero, dropout_nonzero, zeros_proportion * 100),
         xlab = "Expression", breaks = 30, xlim = c(0, max(sim_initial$data)))
  }
  
  par(mfrow = c(1, 1))
}

# Запуск визуализации
cat("\n=== Визуализация эффекта dropout ===\n")
visualize_dropout_effect()