library('GEOquery')
library(biomaRt)
library(ggplot2)
library(hgu133plus2.db)
library(ComplexHeatmap)
library(progress)
library(reshape)
library(plotly)
library(ggrastr)
library(ggpubr)
library(linseed)
library(svglite)
library(DualSimplex)
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
    # Радиус = медиана расстояний до k-го соседа
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
                                     keep_lower = FALSE,
                                     plot = TRUE,          # новый параметр
                                     verbose = TRUE) {
  
  iter <- 1
  repeat {
    # 1. Обновить аннотацию плотности на текущих данных
    add_density_annov2(k = k, genes = genes)
    
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
    
    # 3. Применить фильтр
    filtered_data <- threshold_filter(
      eset = current_data,
      feature = "density",
      threshold = threshold,
      genes = genes,
      keep_lower = keep_lower   
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
n_ct <- 4

library('GEOquery')
dataset <- "GSE11058"
gse <- getGEO(dataset, AnnotGPL = T)
gse <- gse[[1]]
data_raw <- Biobase::exprs(gse)
#linearize the data!
data_raw <- linearize_dataset(data_raw)

# probes which are present in data
probe_mapping <-
  biomaRt::select(hgu133plus2.db::hgu133plus2.db,
                  rownames(data_raw),
                  c("SYMBOL"),
                  keytype = "PROBEID")
probe_mapping <-
  probe_mapping[probe_mapping$PROBEID %in% rownames(data_raw),]
# not empty result gene id (now we  many-to-many in in genes to probes)
probe_mapping <- probe_mapping[!is.na(probe_mapping$SYMBOL),]
## want to collapse gene for each gene and take mean value of data
data_raw <- data_raw[probe_mapping$PROBEID, ]
rownames(data_raw) <- probe_mapping$SYMBOL
# should be unique gene name as row name after this
data_raw <-
  tapply(data_raw, list(row.names(data_raw)[row(data_raw)], colnames(data_raw)[col(data_raw)]), FUN = median)

### True H
pdata <- pData(gse)
component_names <- c("Jurkat", "IM-9", "Raji", "THP-1")

mix_a <- c(2.5, 1.25, 2.5, 3.75)
mix_a <- mix_a / sum(mix_a)
mix_b <- c(0.5, 3.17, 4.75, 1.58)
mix_b <- mix_b / sum(mix_b)
mix_c <- c(0.1, 4.95, 1.65, 3.3)
mix_c <- mix_c/ sum(mix_c)
mix_d <- c(0.02, 3.33, 3.33, 3.33)
mix_d <- mix_d / sum(mix_d)

proportions <- matrix(0 , nrow = dim(data_raw)[[2]], ncol = n_ct )

proportions[1:3, 1] <- 1
proportions[4:6, 2] <- 1
proportions[7:9, 3] <- 1
proportions[10:12, 4] <- 1

proportions[13, ] <- mix_a
proportions[14, ] <- mix_a
proportions[15, ] <- mix_a

proportions[16, ] <- mix_b
proportions[17, ] <- mix_b
proportions[18, ] <- mix_b

proportions[19, ] <- mix_c
proportions[20, ] <- mix_c
proportions[21, ] <- mix_c

proportions[22, ] <- mix_d
proportions[23, ] <- mix_d
proportions[24, ] <- mix_d

rownames(proportions) <- rownames(pdata)
colnames(proportions) <- component_names
true_proportions <- t(proportions)

### True W
components <- lapply(c(1:n_ct), function(comp_num){
  component_subset_columns <- colnames(true_proportions[,true_proportions[comp_num, ] == 1])
  component_subset <- data_raw[, component_subset_columns]
  component_vector <- as.matrix(rowMedians(component_subset))
  rownames(component_vector) <- rownames(component_subset)
  return(component_vector)
  
})
true_basis <- do.call(cbind, components)
colnames(true_basis) <- component_names
true_marker_list <- get_signature_markers(true_basis, n_marker_genes = 100)

print(paste("Dim for the whole data:", toString(dim(data_raw))))
print(paste("Dim for hidden proportions (H):", toString(dim(true_proportions))))
print(paste("Dim for hidden basis (W):", toString(dim(true_basis))))

dso <- DualSimplexSolver$new()
dso$set_data(data_raw)
# remove "RPLS", "LOC", "ORF", "SNOR", keep coding only
dso$basic_filter()
dso$project(n_ct)
print(paste("Dim for the prefiltered data:", toString(dim(dso$get_data()))))

dso$plot_svd(dims = 1:15)
dso$plot_projected("zero_distance", "zero_distance", use_dims = list(3:4))


add_density_annov2(k = 20, genes = T)

density_values <- dso[["st"]][["data"]]@featureData@data[["density"]]
df <- data.frame(density = density_values)

p <- ggplot(df, aes(x = density)) +
  geom_histogram(fill = "lightblue", color = "white", bins = 200) +
  labs(title = "Gene density distribution",
       x = "Density",
       y = "Freq.") +
  theme_minimal()

# Добавить линию порога, если нужно
threshold <- 1
p <- p + geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1)

print(p)

dso <- iterative_density_filter(
  dso, 
  threshold = 1,
  n_cell_types = 4,
  genes = T, 
  k=20, 
  max_iter = 100,
  verbose = TRUE
)

dso$project(4)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(3:4)
)

set.seed(33)
dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()

ptb <- coerce_pred_true_basis(solution$W, true_basis[rownames(solution$W), ])
ptp <- coerce_pred_true_props(solution$H, true_proportions)
plot_ptp_scatter(ptp)
plot_ptb_scatter(ptb)

ptp_lines <- linseed::plotProportions(
  as.data.frame(ptp[[1]]),
  as.data.frame(ptp[[2]]),
  pnames = c("predicted", "true"),
  point_size = 1,
  line_size = 0.7
) + theme_bw(base_size = 12) + theme(
  legend.title = element_blank(),
  legend.position = "bottom",
  axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1)
)

ptp_lines