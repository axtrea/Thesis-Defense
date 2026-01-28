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

# Модифицированная функция для оценки качества с гауссовым шумом
evaluate_gaussian_noise_effect <- function(deviation_vals, additive_vals, n_runs = 3) {
  results <- list()
  all_ptp_plots <- list()  # Список для хранения данных для графиков ptp
  
  for (deviation in deviation_vals) {
    for (additive in additive_vals) {
      cat(sprintf("\n=== Testing: deviation = %.1f, additive = %s ===\n", 
                  deviation, ifelse(additive, "TRUE", "FALSE")))
      
      # Повторяем несколько раз для статистики
      run_metrics <- lapply(1:n_runs, function(run) {
        cat(sprintf("Run %d/%d... ", run, n_runs))
        
        # Создаем симуляцию
        set.seed(run)
        sim_initial <- create_simulation(n_genes = 10000, n_samples = 100, n_cell_types = 3)
        
        # Добавляем гауссов шум
        sim_noisy <- add_gaussian_noise(
          sim_initial, 
          noise_deviation = deviation, 
          additive = additive,
          protect_genes = c(),  # Можно добавить защиту определенных генов
          protect_samples = c() # Можно добавить защиту определенных образцов
        )
        
        # Проверяем нулевые строки (опционально)
        zero_rows_final <- sum(rowSums(sim_noisy$data) == 0)
        if (zero_rows_final > 0) {
          cat(sprintf("Warning: %d zero rows after noise addition!\n", zero_rows_final))
        }
        
        tryCatch({
          # Запускаем DualSimplex
          dso <- DualSimplexSolver$new()
          dso$set_data(sim_noisy$data, max_dim = 10)
          dso$project(3)
          dso$init_solution("random")
          dso$default_optimization()
          solution <- dso$finalize_solution()
          result_W <- solution$W
          result_H <- solution$H
          true_H = sim_initial$proportions
          true_basis = sim_initial$basis
          
          # Нормализация H
          norm_res_list <- lapply(1:dim(result_H)[[2]], function(col_ind) {
            cur_col <- result_H[, col_ind]
            normalized_column <- (cur_col - min(cur_col)) / (max(cur_col) - min(cur_col))
            normalized_column <- cur_col / sum(cur_col)
            return(normalized_column)
          })
          normalized_H <- do.call(cbind, norm_res_list)
          
          # Получаем данные для ptp графика
          ptp <- coerce_pred_true_props(normalized_H, true_H)
          
          # Сохраняем данные для этого запуска
          ptp_data <- list(
            pred_H = ptp[[1]],
            true_H = ptp[[2]],
            noise_deviation = deviation,
            additive = additive,
            run = run
          )
          
          # Добавляем в список всех графиков
          plot_id <- paste0("dev", deviation, "_add", as.numeric(additive), "_run", run)
          all_ptp_plots[[plot_id]] <<- ptp_data
          
          # Сравниваем с истинными пропорциями
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
            noise_deviation = deviation,
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
            noise_deviation = deviation,
            additive = additive,
            run = run,
            pearson_mean = NA,
            rmse_mean = NA
          ))
        })
      })
      
      results[[paste0("dev", deviation, "_add", as.numeric(additive))]] <- do.call(rbind, run_metrics)
    }
  }
  
  # Объединяем все результаты
  final_results <- do.call(rbind, results)
  rownames(final_results) <- NULL
  
  # Возвращаем и результаты метрик, и данные для графиков
  return(list(
    metrics = final_results,
    ptp_plots_data = all_ptp_plots
  ))
}

# Запускаем оценку с разными параметрами гауссова шума
deviation_vals <- c(0.1, 0.5, 1, 2, 5)  # Стандартные отклонения
additive_vals <- c(TRUE, FALSE)         # Аддитивный и мультипликативный шум

evaluation_results <- evaluate_gaussian_noise_effect(deviation_vals, additive_vals, n_runs = 1)
results <- evaluation_results$metrics
all_ptp_plots <- evaluation_results$ptp_plots_data

# Убираем NA значения для визуализации
results_clean <- na.omit(results)

# Функция для создания комбинированного графика
create_combined_ptp_plot <- function(ptp_plots_data) {
  if (length(ptp_plots_data) == 0) {
    cat("Нет данных для построения графиков.\n")
    return(NULL)
  }
  
  # Создаем список графиков
  plot_list <- list()
  
  for (i in seq_along(ptp_plots_data)) {
    plot_data <- ptp_plots_data[[i]]
    plot_name <- names(ptp_plots_data)[i]
    
    # Создаем данные для scatter plot
    plot_df <- data.frame()
    for (cell_type in 1:nrow(plot_data$pred_H)) {
      temp_df <- data.frame(
        Predicted = plot_data$pred_H[cell_type, ],
        True = plot_data$true_H[cell_type, ],
        CellType = paste0("CellType ", cell_type)
      )
      plot_df <- rbind(plot_df, temp_df)
    }
    
    # Создаем scatter plot
    p <- ggplot(plot_df, aes(x = True, y = Predicted, color = CellType)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      facet_wrap(~CellType, scales = "free") +
      labs(
        title = paste0("Noise: σ=", plot_data$noise_deviation, 
                       ", ", ifelse(plot_data$additive, "Additive", "Multiplicative"),
                       " (Run ", plot_data$run, ")"),
        x = "True Proportions",
        y = "Predicted Proportions"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)
      )
    
    plot_list[[plot_name]] <- p
  }
  
  # Комбинируем все графики в одну фигуру
  n_plots <- length(plot_list)
  n_cols <- 2  # Количество колонок в сетке
  n_rows <- ceiling(n_plots / n_cols)
  
  # Используем ggarrange для компоновки
  combined_plot <- ggarrange(
    plotlist = plot_list,
    ncol = n_cols,
    nrow = n_rows,
    common.legend = TRUE,
    legend = "bottom"
  )
  
  # Добавляем общий заголовок
  combined_plot <- annotate_figure(
    combined_plot,
    top = text_grob("Predicted vs True Proportions for Different Noise Levels", 
                    face = "bold", size = 14)
  )
  
  return(combined_plot)
}

# Создаем комбинированный график
if (length(all_ptp_plots) > 0) {
  combined_ptp_plot <- create_combined_ptp_plot(all_ptp_plots)
  print(combined_ptp_plot)
  
  # Можно сохранить график
  # ggsave("combined_ptp_plots.png", combined_ptp_plot, width = 16, height = 20, dpi = 300)
  # ggsave("combined_ptp_plots.pdf", combined_ptp_plot, width = 16, height = 20)
}

# Остальной код для визуализации метрик остается без изменений
if (nrow(results_clean) > 0) {
  # Преобразуем additive в фактор для лучшей визуализации
  results_clean$additive <- factor(results_clean$additive, 
                                   levels = c(TRUE, FALSE),
                                   labels = c("Additive", "Multiplicative"))
  
  # График для Pearson correlation
  p1 <- ggplot(results_clean, aes(x = factor(noise_deviation), y = pearson_mean, 
                                  color = additive, group = additive)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Noise Deviation (σ)", y = "Pearson Correlation", 
         color = "Noise Type", 
         title = "Влияние гауссова шума на точность деконволюции (Pearson)") +
    theme_minimal() +
    ylim(-0.1, 1.1) +
    scale_color_manual(values = c("Additive" = "blue", "Multiplicative" = "red"))
  
  # График для RMSE
  p2 <- ggplot(results_clean, aes(x = factor(noise_deviation), y = rmse_mean,
                                  color = additive, group = additive)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    labs(x = "Noise Deviation (σ)", y = "RMSE", 
         color = "Noise Type",
         title = "Влияние гауссова шума на точность деконволюции (RMSE)") +
    theme_minimal() +
    scale_color_manual(values = c("Additive" = "blue", "Multiplicative" = "red"))
  
  # Показываем графики
  print(p1)
  print(p2)
  
  # ... остальной код визуализации метрик ...
} else {
  cat("\nВсе запуски завершились с ошибками. Попробуйте уменьшить параметры шума.\n")
}
  
  
  # Сводная статистика
  cat("\n=== Сводная статистика ===\n")
  print(agg_data)
  
  # Дополнительный график: сравнение двух типов шума
  p_compare <- ggplot(results_clean, aes(x = factor(noise_deviation), y = pearson_mean, fill = additive)) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
    labs(x = "Noise Deviation (σ)", y = "Pearson Correlation", 
         fill = "Noise Type",
         title = "Сравнение аддитивного и мультипликативного гауссова шума") +
    theme_minimal() +
    scale_fill_manual(values = c("Additive" = "blue", "Multiplicative" = "red"))
  
  print(p_compare)
  
} else {
  cat("\nВсе запуски завершились с ошибками. Попробуйте уменьшить параметры шума.\n")
}
