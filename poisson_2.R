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

# Модифицированная функция для добавления Пуассоновского шума (для реальных данных)
add_poisson_noise_real <- function(data, scaling_factor = 1,
                                   protect_genes = c(), protect_samples = c()) {
  
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
  
  return(noisy_data)
}

# Функция для сопоставления предсказанных и истинных типов клеток
match_cell_types <- function(pred_H, true_H) {
  n_components <- nrow(pred_H)
  
  # Нормализуем обе матрицы (сумма по столбцам = 1)
  pred_H_norm <- apply(pred_H, 2, function(x) x / sum(x))
  true_H_norm <- apply(true_H, 2, function(x) x / sum(x))
  
  # Сопоставляем типы клеток (строки) используя корреляцию
  correlation_matrix <- matrix(0, n_components, n_components)
  for (i in 1:n_components) {
    for (j in 1:n_components) {
      correlation_matrix[i, j] <- cor(pred_H_norm[i, ], true_H_norm[j, ], 
                                      method = "pearson", use = "complete.obs")
    }
  }
  
  # Находим наилучшее соответствие (максимизируем сумму корреляций)
  matches <- rep(0, n_components)
  used_pred <- rep(FALSE, n_components)
  used_true <- rep(FALSE, n_components)
  
  for (k in 1:n_components) {
    # Находим максимальную корреляцию среди неиспользованных
    max_corr <- -Inf
    best_i <- 0
    best_j <- 0
    
    for (i in 1:n_components) {
      if (!used_pred[i]) {
        for (j in 1:n_components) {
          if (!used_true[j]) {
            if (correlation_matrix[i, j] > max_corr) {
              max_corr <- correlation_matrix[i, j]
              best_i <- i
              best_j <- j
            }
          }
        }
      }
    }
    
    matches[best_i] <- best_j
    used_pred[best_i] <- TRUE
    used_true[best_j] <- TRUE
  }
  
  # Переставляем строки pred_H в соответствии с найденным соответствием
  pred_H_matched <- pred_H_norm[matches, ]
  rownames(pred_H_matched) <- rownames(true_H_norm)
  
  return(pred_H_matched)
}

# Функция для оценки качества на реальных данных
evaluate_noise_effect_real <- function(log_data, linearize_function, true_proportions, 
                                       scaling_vals, 
                                       n_runs = 3, n_components = 3) {
  results <- list()
  
  for (scale in scaling_vals) {
    cat(sprintf("\n=== Testing: scaling_factor = %.1f ===\n", scale))
    
    # Повторяем несколько раз для статистики
    run_metrics <- lapply(1:n_runs, function(run) {
      cat(sprintf("Run %d/%d... ", run, n_runs))
      
      # Добавляем Пуассоновский шум к log-данным
      set.seed(run)
      noisy_log_data <- add_poisson_noise_real(
        log_data, 
        scaling_factor = scale,
        protect_genes = c(),
        protect_samples = c()
      )
      
      # Линеризируем зашумленные данные
      tryCatch({
        noisy_linear_data <- linearize_function(noisy_log_data)
        
        # Запускаем DualSimplex
        dso <- DualSimplexSolver$new()
        dso$set_data(noisy_linear_data, max_dim = 30)
        dso$project(n_components)
        dso$init_solution("random")
        dso$default_optimization()
        solution <- dso$finalize_solution()
        
        # Сопоставляем типы клеток
        pred_H_matched <- match_cell_types(solution$H, true_proportions)
        true_H_norm <- apply(true_proportions, 2, function(x) x / sum(x))
        
        # Вычисляем метрики для каждого типа клеток
        metrics <- sapply(1:n_components, function(i) {
          pred <- pred_H_matched[i, ]
          true <- true_H_norm[i, ]
          c(
            pearson = ifelse(sd(pred, na.rm = TRUE) > 0 & sd(true, na.rm = TRUE) > 0, 
                             cor(pred, true, method = "pearson", use = "complete.obs"), 
                             NA),
            rmse = sqrt(mean((pred - true)^2, na.rm = TRUE))
          )
        })
        
        # Усредняем по типам клеток
        result_df <- data.frame(
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
          scaling_factor = scale,
          run = run,
          pearson_mean = NA,
          rmse_mean = NA
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

# Загрузка и подготовка реальных данных
library('GEOquery')
library(biomaRt)

# Определяем функцию для линеаризации (если не существует)
if (!exists("linearize_dataset")) {
  linearize_dataset <- function(log_data) {
    # Экспоненциальное преобразование для возврата к линейной шкале
    # с учетом что входные данные в log2 шкале (log2(x+1))
    linear_data <- 2^log_data - 1
    # Защита от отрицательных значений из-за численных погрешностей
    linear_data[linear_data < 0] <- 0
    return(linear_data)
  }
}

n_ct <- 3
dataset <- "GSE19830"
gse <- getGEO(dataset, AnnotGPL = T)
gse <- gse[[1]]
data_raw <- Biobase::exprs(gse)

# Преобразуем в log2 шкалу (если еще не в log)
# Предполагаем, что данные могут быть в линейной шкале
# Проверяем по максимальному значению
max_val <- max(data_raw)
if (max_val > 1000) {
  # Данные в линейной шкале, преобразуем в log2
  data_log <- log2(data_raw + 1)
} else {
  # Данные уже в log шкале
  data_log <- data_raw
}

# Маппинг probe IDs на gene symbols
probe_mapping <- biomaRt::select(rat2302.db::rat2302.db, rownames(data_log), c("SYMBOL", "GENETYPE"))
probe_mapping <- probe_mapping[probe_mapping$PROBEID %in% rownames(data_log),]
probe_mapping <- probe_mapping[!is.na(probe_mapping$SYMBOL),]
probe_mapping <- probe_mapping[!duplicated(probe_mapping$PROBEID),]

data_log <- data_log[probe_mapping$PROBEID, ]
rownames(data_log) <- probe_mapping$SYMBOL
data_log <- tapply(data_log, list(row.names(data_log)[row(data_log)], 
                                  colnames(data_log)[col(data_log)]), FUN = median)

# Линеаризуем данные для вычисления истинных пропорций
data_linear <- linearize_dataset(data_log)

# Истинные пропорции (H) - вычисляем на линейных данных
component_names <- c("Liver", "Brain", "Lung")
pdata <- pData(gse)
parsed_proportions <- strsplit(pdata$source_name_ch1,split = '/')
parsed_proportions <- lapply(parsed_proportions, function(sample_props){
  return(strtoi(gsub("\\D*(\\d+)\\D*","\\1",sample_props)))
})
true_proportions <- matrix(unlist(parsed_proportions), ncol=3, byrow=TRUE)/100
rownames(true_proportions) <- rownames(pdata)
colnames(true_proportions) <- component_names
true_proportions <- t(true_proportions)

# Истинные сигнатуры (W) на линейных данных
components <- lapply(1:n_ct, function(comp_num){
  component_subset_columns <- colnames(true_proportions[,true_proportions[comp_num, ] == 1])
  component_subset <- data_linear[, component_subset_columns]
  component_vector <- as.matrix(rowMedians(component_subset))
  rownames(component_vector) <- rownames(component_subset)
  return(component_vector)
})
true_basis <- do.call(cbind, components)
colnames(true_basis) <- component_names

# Используем только coding genes для ускорения вычислений
coding_gene_info = probe_mapping[probe_mapping$GENETYPE == "protein-coding",]
gene_anno <- list(CODING = coding_gene_info$SYMBOL)

# Берем подмножество данных для ускорения
coding_genes <- rownames(data_log)[rownames(data_log) %in% gene_anno$CODING]
data_log_coding <- data_log[coding_genes, ]

# Запускаем оценку с разными параметрами шума
scaling_vals <- c(0.1, 0.5, 1, 2, 5)

results <- evaluate_noise_effect_real(
  log_data = data_log_coding,
  linearize_function = linearize_dataset,
  true_proportions = true_proportions,
  scaling_vals = scaling_vals,
  n_runs = 1,
  n_components = n_ct
)

# Убираем NA значения для визуализации
results_clean <- na.omit(results)

# Визуализация результатов
if (nrow(results_clean) > 0) {
  # График для Pearson correlation
  p1 <- ggplot(results_clean, aes(x = factor(scaling_factor), y = pearson_mean)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1, group = 1) +
    labs(x = "Scaling Factor", y = "Pearson Correlation", 
         title = "Влияние Пуассоновского шума на точность деконволюции (Pearson)") +
    theme_minimal() +
    ylim(-0.1, 1.1)
  
  # График для RMSE
  p2 <- ggplot(results_clean, aes(x = factor(scaling_factor), y = rmse_mean)) +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1, group = 1) +
    labs(x = "Scaling Factor", y = "RMSE",
         title = "Влияние Пуассоновского шума на точность деконволюции (RMSE)") +
    theme_minimal()
  
  # Показываем графики
  print(p1)
  print(p2)
  
  # Сводная статистика
  cat("\n=== Сводная статистика ===\n")
  agg_data <- aggregate(cbind(pearson_mean, rmse_mean) ~ scaling_factor, 
                        data = results_clean, FUN = mean, na.action = na.pass)
  print(agg_data)
  
  # График для сравнения метрик
  library(reshape2)
  results_melted <- reshape2::melt(results_clean, id.vars = c("scaling_factor", "run"),
                                   measure.vars = c("pearson_mean", "rmse_mean"))
  
  p3 <- ggplot(results_melted, aes(x = factor(scaling_factor), y = value, fill = variable)) +
    geom_boxplot() +
    labs(x = "Scaling Factor", y = "Значение метрики", 
         title = "Распределение метрик качества для разных уровней шума",
         fill = "Метрика") +
    theme_minimal() +
    facet_wrap(~variable, scales = "free_y") +
    theme(legend.position = "none")
  
  print(p3)
  
  # Дополнительная визуализация: сравнение исходных и зашумленных данных
  cat("\n=== Визуализация эффекта шума ===\n")
  
  # Пример для одного значения scaling_factor
  scale_example <- 1
  set.seed(123)
  noisy_log_example <- add_poisson_noise_real(data_log_coding[, 1:20], scaling_factor = scale_example)
  noisy_linear_example <- linearize_dataset(noisy_log_example)
  
  # Создаем data.frame для визуализации
  plot_data <- data.frame(
    Value = c(as.vector(data_linear[coding_genes, 1:20]), 
              as.vector(noisy_linear_example)),
    Type = rep(c("Original", "Noisy"), each = length(as.vector(data_linear[coding_genes, 1:20]))),
    Dataset = rep(rep(1:20, each = length(coding_genes)), 2)
  )
  
  p4 <- ggplot(plot_data, aes(x = Value, fill = Type)) +
    geom_density(alpha = 0.5) +
    scale_x_log10() +
    labs(x = "Expression (log10 scale)", y = "Density",
         title = "Распределение выражений: исходные vs зашумленные данные",
         fill = "Data Type") +
    theme_minimal()
  
  print(p4)
  
} else {
  cat("\nВсе запуски завершились с ошибками. Попробуйте уменьшить параметры шума.\n")
}
