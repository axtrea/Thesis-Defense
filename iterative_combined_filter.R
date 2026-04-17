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

compute_extremeness <- function(coords, k = 10) {
  # Центроид
  centroid <- colMeans(coords)
  
  # Расстояние до центроида
  dist_to_centroid <- sqrt(rowSums((coords - matrix(centroid, 
                                                    nrow(coords), 
                                                    ncol(coords), 
                                                    byrow = TRUE))^2))
  
  # kNN расстояния
  knn <- dbscan::kNN(coords, k = k)
  knn_dist <- rowMeans(knn$dist, na.rm = TRUE)
  
  # Нормализация
  ext1 <- scale(dist_to_centroid)[,1]
  ext2 <- scale(knn_dist)[,1]
  
  extremeness <- ext1 + ext2
  
  names(extremeness) <- rownames(coords)
  return(extremeness)
}

add_extremeness_anno <- function(dso, k = 10, genes = TRUE) {
  
  eset <- dso$get_data()
  
  # координаты
  coords <- if (genes) dso$st$proj$X else dso$st$proj$Omega
  
  extremeness <- compute_extremeness(coords, k = k)
  
  anno <- get_anno(eset, genes)
  anno$extremeness <- extremeness[rownames(anno)]
  
  eset <- set_anno(anno, eset, genes)
  
  dso$set_data(eset)
  
  return(dso)
}

compute_extremeness <- function(coords, k = 10) {
  # Центроид
  centroid <- colMeans(coords)
  
  # Расстояние до центроида
  dist_to_centroid <- sqrt(rowSums((coords - matrix(centroid, 
                                                    nrow(coords), 
                                                    ncol(coords), 
                                                    byrow = TRUE))^2))
  
  # kNN расстояния
  knn <- dbscan::kNN(coords, k = k)
  knn_dist <- rowMeans(knn$dist, na.rm = TRUE)
  
  # Нормализация
  ext1 <- scale(dist_to_centroid)[,1]
  ext2 <- scale(knn_dist)[,1]
  
  extremeness <- ext1 + ext2
  
  names(extremeness) <- rownames(coords)
  return(extremeness)
}
#' Комбинированная итеративная фильтрация по плотности и экстремальности
#'
#' @param dso Объект данных (класс DualSimplexSolver)
#' @param density_threshold Порог плотности: удаляются объекты с density < density_threshold
#' @param extremeness_threshold Порог экстремальности: удаляются объекты с extremeness < extremeness_threshold
#' @param k Количество соседей для расчёта экстремальности (передаётся в add_extremeness_anno)
#' @param genes Логический: TRUE – фильтровать гены, FALSE – образцы
#' @param max_iter Максимальное число итераций
#' @param verbose Выводить ли сообщения о прогрессе
#'
#' @return Обновлённый объект dso (невидимо)
iterative_combined_filter <- function(dso,
                                      density_threshold = 1,
                                      extremeness_threshold = 0,
                                      k = 10,
                                      genes = FALSE,
                                      n_cell_types = 4,
                                      max_iter = 100,
                                      verbose = TRUE) {
  
  iter <- 1
  repeat {
    # 1. Добавить аннотации плотности и экстремальности к текущим данным
    dso$add_density_anno(genes = genes)          
    add_extremeness_anno(dso, k = k, genes = genes)  
    
    # 2. Получить текущие данные
    current_data <- dso$get_data()
    
    # 3. Применить фильтрацию по обоим условиям
    anno <- get_anno(current_data, genes)        # извлекаем аннотацию
    
    # Логический вектор: оставляем только те объекты, у которых оба показателя >= порогов
    keep <- anno$density >= density_threshold & anno$extremeness >= extremeness_threshold
    
    # Отфильтровать данные
    if (genes) {
      filtered_data <- current_data[keep, , drop = FALSE]
    } else {
      filtered_data <- current_data[, keep, drop = FALSE]
    }
    
    # 4. Подсчитать изменения
    n_before <- if (genes) nrow(current_data) else ncol(current_data)
    n_after  <- if (genes) nrow(filtered_data) else ncol(filtered_data)
    removed <- n_before - n_after
    
    if (verbose) {
      cat(sprintf("Итерация %d: удалено %d объектов (density < %g или extremeness < %g)\n", 
                  iter, removed, density_threshold, extremeness_threshold))
    }
    
    # 5. Условия остановки
    if (removed == 0) {
      if (verbose) cat("Объектов, не проходящих фильтр, не осталось. Остановка.\n")
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
    
    # 6. Обновить данные и перепроецировать
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

#dso$add_density_anno(genes = TRUE)          
#add_extremeness_anno(dso, k = 10, genes = TRUE)
#current_data <- dso$get_data()
#anno <- get_anno(current_data, genes = TRUE)
#keep <- anno$density >= 100 & anno$extremeness >= 0


dso <- iterative_combined_filter(dso,
                                 density_threshold = 1,
                                 extremeness_threshold = -0.5,
                                 k = 10,
                                 genes = TRUE,
                                 max_iter = 1,
                                 verbose = TRUE)

dso$project(4)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)



set.seed(1)
dso$init_solution("random")
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  with_solution = TRUE,
  use_dims = list(2:3)
)
dso$default_optimization()
solution <- dso$finalize_solution()

cur_H <- solution$H

norm_res_list <- lapply(1:dim(cur_H)[[2]], function(col_ind) {
  cur_col <- cur_H[, col_ind]
  normalized_column <-
    (cur_col - min(cur_col)) / (max(cur_col) - min(cur_col))
  normalized_column <- cur_col / sum(cur_col)
  return(normalized_column)
})
normalized_H <- do.call(cbind, norm_res_list)

ptp <- coerce_pred_true_props(normalized_H, true_proportions)
plot_ptp_scatter(ptp)
ptb <-
  coerce_pred_true_basis(solution$W, true_basis[rownames(solution$W),])
plot_ptb_scatter(ptb)
#Visualize as lines
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