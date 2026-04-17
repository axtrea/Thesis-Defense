library("DualSimplex")
library(dplyr)
library(reshape2) 
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(RColorBrewer)

#' Итеративная фильтрация по плотности
#'
#' @param dso Объект данных 
#' @param threshold Порог плотности (удаляются объекты с density < threshold)
#' @param genes Логический: TRUE для генов, FALSE для образцов (по умолчанию образцы)
#' @param density_radius Радиус для расчёта плотности (передаётся в add_density_anno)
#' @param max_iter Максимальное число итераций
#' @param verbose Выводить ли сообщения о прогрессе
#'
#' @return Обновлённый объект dso (невидимо)
iterative_density_filter <- function(dso, 
                                     threshold = 1,
                                     n_cell_types = NULL,
                                     genes = FALSE, 
                                     density_radius = NULL,
                                     plot = TRUE,
                                     max_iter = 100,
                                     verbose = TRUE) {
  
  
  iter <- 1
  repeat {
    # 1. Обновить аннотацию плотности на текущих данных
    dso$add_density_anno(radius = density_radius, genes = genes)
    
    # 2. Получить текущие данные
    current_data <- dso$get_data()

    # --- Построение гистограммы, если нужно ---
    if (plot) {
      # Извлекаем аннотацию (featureData или phenoData)
      anno <- get_anno(current_data, genes)
      density_values <- anno$density
      
      # Создаём датафрейм для ggplot
      df <- data.frame(density = density_values)
      
      # Строим гистограмму
      p <- ggplot(df, aes(x = density)) +
        geom_histogram(fill = "lightblue", color = "white", bins = 100) +
        geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
        labs(title = paste("Gene density distribution - iteration", iter),
             x = "Density", y = "Freq.") +
        theme_minimal()
      
      print(p)   # выводим на текущее графическое устройство
    }
    
    # 3. Применить фильтр: оставить объекты с density >= threshold
    filtered_data <- threshold_filter(
      eset = current_data,
      feature = "density",
      threshold = threshold,
      genes = genes,
      keep_lower = FALSE   # сохраняем те, у которых значение >= threshold
    )
    
    # 4. Подсчитать изменения
    n_before <- if (genes) nrow(current_data) else ncol(current_data)
    n_after  <- if (genes) nrow(filtered_data) else ncol(filtered_data)
    removed <- n_before - n_after
    
    if (verbose) {
      cat(sprintf("Итерация %d: удалено %d объектов (density < %g)\n", 
                  iter, removed, threshold))
    }
    
    # 5. Условия остановки
    if (removed == 0) {
      if (verbose) cat("Объектов ниже порога не осталось. Остановка.\n")
      break
    }
    if (n_after == 0) {
      warning("Все объекты удалены. Остановка.")
      dso$set_data(filtered_data)
      break
    }
    if (iter >= max_iter) {
      warning("Достигнуто максимальное число итераций без сходимости.")
      dso$set_data(filtered_data)
      break
    }
    
    # 6. Обновить данные в dso для следующей итерации
    dso$set_data(filtered_data)
    dso$project(n_cell_types)
    
    iter <- iter + 1
  }
  
  invisible(dso)
}


#Simulation
n_ct = 3
set.seed(3)
sim <- create_simulation(n_genes = 10000,
                         n_samples = 100,
                         n_cell_types = n_ct,
                         with_marker_genes = FALSE)

true_basis <- sim$basis
true_proportions <- sim$proportions


sim <- sim %>% add_noise(noise_deviation = 5)
data_raw <- sim$data
dso <- DualSimplexSolver$new()
dso$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
dso$project(n_ct)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

dso$add_density_anno(genes = T)
density_values <- dso[["st"]][["data"]]@featureData@data[["density"]]
df <- data.frame(density = density_values)

p <- ggplot(df, aes(x = density)) +
  geom_histogram(fill = "lightblue", color = "white", bins = 100) +
  labs(title = "Gene density distribution",
       x = "Density",
       y = "Freq.") +
  theme_minimal()

# Добавить линию порога, если нужно
threshold <- 1
p <- p + geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1)

print(p)

#Filter
dso <- iterative_density_filter(
  dso, 
  threshold = 100,
  n_cell_types = n_ct,
  genes = T, 
  density_radius = NULL, 
  max_iter = 100,
  verbose = TRUE
)

dso$project(n_ct)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

#Solution
dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()
ptp <- coerce_pred_true_props(solution$H, true_proportions)
plot_ptp_scatter(ptp)
plot_ptp_lines(ptp)
ptb <- coerce_pred_true_basis(solution$W, true_basis[rownames(solution$W),])
plot_ptb_scatter(ptb)


#History of solution
colnames(dso$st$proj$X) <- c("R1", "R2", "R3")
colnames(dso$st$proj$Omega) <- c("S1", "S2", "S3")
proj_solution_history <-
  dso$plot_projected(use_dims = (2:3),
                     wrap = F,
                     with_legend = F)
proj_solution_history[[1]] # genes in a space of samples
proj_solution_history[[2]] # samples in a space of genes

common_samples <- intersect(colnames(solution$H), colnames(true_proportions))

# Фильтруем true_proportions, оставляя только общие образцы
true_proportions_filtered <- true_proportions[, common_samples, drop = FALSE]

# Проверяем размерность
dim(true_proportions_filtered)
ptp <- coerce_pred_true_props(solution$H, true_proportions_filtered)
plot_ptp_scatter(ptp)
plot_ptp_lines(ptp)