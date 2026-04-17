library("DualSimplex")
library(dplyr)
library(reshape2) 
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(RColorBrewer)
library(dbscan)

add_density_annov2= function(k = 10, genes = T) {
  dso$st$data <- add_density_annotationv2(
    dso$st$data,
    dso$st$proj,
    k = k,
    genes = genes
  )
}
add_density_annotationv2 <- function(eset, proj, genes = TRUE, radius = NULL, k = 10) {
  # Получаем текущую аннотацию (featureData или phenoData)
  anno <- get_anno(eset, genes)
  
  # Выбираем координаты проекции
  if (genes) {
    coords <- proj$X
  } else {
    coords <- proj$Omega
  }
  
  # Если радиус не задан, вычисляем его по k-дистанциям
  if (is.null(radius)) {
    # Находим k ближайших соседей для каждой точки
    knn_result <- dbscan::kNN(coords, k = k)
    # Радиус = медиана расстояний до k-го соседа (можно заменить на любой квантиль)
    radius <- median(knn_result$dist[, k], na.rm = TRUE)
    # Альтернативно: radius <- quantile(knn_result$dist[, k], 0.95, na.rm = TRUE)
  }
  
  # Находим всех соседей в пределах вычисленного радиуса
  nn_result <- dbscan::frNN(coords, eps = radius)
  
  # Количество соседей (включая саму точку)
  nn_count <- unlist(lapply(nn_result$id, length))
  # Среднее расстояние до соседей
  nn_mean_distance <- unlist(lapply(nn_result$dist, mean))
  
  # Присваиваем имена строк для безопасного сопоставления
  names(nn_count) <- rownames(coords)
  names(nn_mean_distance) <- rownames(coords)
  
  # Обновляем аннотацию
  anno$density <- nn_count[rownames(anno)]
  anno$mean_nn_distance <- nn_mean_distance[rownames(anno)]
  
  eset <- set_anno(anno, eset, genes)
  return(eset)
}
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
                                     k = 10,
                                     max_iter = 100,
                                     verbose = TRUE) {
  
  
  iter <- 1
  repeat {
    # 1. Обновить аннотацию плотности на текущих данных
    add_density_annov2( k = k, genes = genes)
    
    # 2. Получить текущие данные
    current_data <- dso$get_data()
    
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
set.seed(3)
sim <- create_simulation(n_genes = 10000,
                         n_samples = 100,
                         n_cell_types = 3,
                         with_marker_genes = FALSE)

true_basis <- sim$basis
true_proportions <- sim$proportions


sim <- sim %>% add_noise(noise_deviation = 2)
data_raw <- sim$data
dso <- DualSimplexSolver$new()
dso$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
dso$project(3)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

dso$plot_projected(
  "plane_distance",
  "plane_distance",
  use_dims = list(2:3)
)

add_density_annov2(k = 1000, genes = T)
density_values <- dso[["st"]][["data"]]@featureData@data[["density"]]
df <- data.frame(density = density_values)

p <- ggplot(df, aes(x = density)) +
  geom_histogram(fill = "lightblue", color = "white", bins = 100) +
  labs(title = "Gene density distribution",
       x = "Density",
       y = "Freq.") +
  theme_minimal()

# Добавить линию порога, если нужно
threshold <- 10
p <- p + geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1)

print(p)

#Filter
dso <- iterative_density_filter(
  dso, 
  threshold = 10,
  n_cell_types = 3,
  genes = T, 
  k=1000, 
  max_iter = 100,
  verbose = TRUE
)

dso$project(3)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

dso$plot_projected(
  "plane_distance",
  "plane_distance",
  use_dims = list(2:3)
)


colors <- which_marker(rownames(fData(dso$st$data)), dso$st$marker_genes)
plot_markers <-
  plot_projection_points(
    dso$st$proj,
    use_dims = (2:3),
    spaces = c("X"),
    pt_size = 1,
    color = colors
  ) +
  scale_color_manual(values = colors_v[3:5],
                     na.value = adjustcolor("grey70", alpha.f = 0.7)) +
  labs(col = "Marker Cell Type", x = "R2", y = "R3") +
  theme_bw(base_family = "sans", base_size = 12) +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )


plot_markers



#Solution
dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()
ptp <- coerce_pred_true_props(solution$H, true_proportions)
plot_ptp_scatter(ptp)
plot_ptp_lines(ptp)
ptb <- coerce_pred_true_basis(solution$W, true_basis[rownames(solution$W),])
plot_ptb_scatter(ptb)