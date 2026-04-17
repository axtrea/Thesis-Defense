library(biomaRt)
library('GEOquery')
library(ggplot2)
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
library(hgu133plus2.db)

# ------------------------------------------------------------
# Функция-обёртка для добавления плотности (исправлена)
# ------------------------------------------------------------
add_annov2 <- function(n_ct = 3, genes = TRUE) {
  dso$st$data <- add_annotationv2(
    eset = dso$st$data,
    proj = dso$st$proj,
    genes = genes,
    n_ct = n_ct
  )
}

# ------------------------------------------------------------
# Основная функция расчёта плотности с фиксированным радиусом на основе MAD
# ------------------------------------------------------------
add_annotationv2 <- function(eset, proj, genes = TRUE, radius = NULL, n_ct = 3) {
  anno <- get_anno(eset, genes)
  
  # Выбираем координаты проекции
  if (genes) {
    coords <- proj$X
  } else {
    coords <- proj$Omega
  }
  
  # Если радиус не задан, вычисляем его как MAD всех координат (кроме первой) * sqrt(n_ct)
  if (is.null(radius)) {
    # coords может быть матрицей или датафреймом, берем все столбцы (если coords одномерный, то берем его)
    if (is.matrix(coords) || is.data.frame(coords)) {
      if (ncol(coords) > 1) {
        coords_mad <- coords[, 2:ncol(coords), drop = FALSE]
      } else {
        coords_mad <- coords
      }
      radius <- stats::mad(coords_mad) * sqrt(n_ct)
    } else {
      radius <- stats::mad(coords) * sqrt(n_ct)
    }
  }
  
  # Находим соседей в радиусе
  nn_result <- dbscan::frNN(coords, eps = radius)
  nn_count <- unlist(lapply(nn_result$id, length))
  nn_mean_distance <- unlist(lapply(nn_result$dist, mean))
  
  names(nn_count) <- rownames(coords)
  names(nn_mean_distance) <- rownames(coords)
  
  anno$density <- nn_count[rownames(anno)]
  anno$mean_nn_distance <- nn_mean_distance[rownames(anno)]
  
  eset <- set_anno(anno, eset, genes)
  return(eset)
}

# ------------------------------------------------------------
# Итеративная фильтрация по плотности
# ------------------------------------------------------------
iterative_density_filter <- function(dso, 
                                     threshold = 1,
                                     n_cell_types = 3,   # задаём значение по умолчанию
                                     genes = FALSE, 
                                     max_iter = 100,
                                     verbose = TRUE) {
  
  iter <- 1
  repeat {
    # 1. Обновить аннотацию плотности на текущих данных
    add_annov2(n_ct = n_cell_types, genes = genes)
    
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

# ------------------------------------------------------------
# Пример использования (Simulation)
# ------------------------------------------------------------
set.seed(3)
n_ct <- 15
sim <- create_simulation(n_genes = 10000,
                         n_samples = 100,
                         n_cell_types = n_ct,
                         with_marker_genes = FALSE)

true_basis <- sim$basis
true_proportions <- sim$proportions

sim <- sim %>% add_noise(noise_deviation = 3)
data_raw <- sim$data
dso <- DualSimplexSolver$new()
dso$set_data(data_raw, max_dim = 30, sinkhorn_tol = 1e-17)
dso$project(n_ct)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

# Добавляем аннотацию плотности для визуализации (до фильтрации)
add_annov2(n_ct = n_ct, genes = TRUE)

# Извлекаем значения плотности (не extremeness, так как её нет в этой версии)
# Для гистограммы используем density
density_values <- dso[["st"]][["data"]]@featureData@data[["density"]]
df <- data.frame(density = density_values)

p <- ggplot(df, aes(x = density)) +
  geom_histogram(fill = "lightblue", color = "white", bins = 200) +
  labs(title = "Gene density distribution",
       x = "Density",
       y = "Freq.") +
  theme_minimal() +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", size = 1)
print(p)

# Запуск итеративной фильтрации
dso <- iterative_density_filter(dso, 
                                threshold = 1, 
                                n_cell_types = n_ct, 
                                genes = TRUE)

# Повторная проекция и визуализация
dso$project(n_ct)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)

# Оптимизация
dso$init_solution("random")
dso$default_optimization()
solution <- dso$finalize_solution()
ptp <- coerce_pred_true_props(solution$H, true_proportions)
plot_ptp_lines(ptp)





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

dso$plot_projected("zero_distance", "zero_distance", use_dims = list(3:4))


# Запуск итеративной фильтрации
dso <- iterative_density_filter(dso, 
                                threshold = 1, 
                                n_cell_types = n_ct, 
                                genes = TRUE)

# Повторная проекция и визуализация
dso$project(n_ct)
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  use_dims = list(2:3)
)



#set.seed(1)
dso$init_solution("random")
dso$default_optimization()


#dso$init_solution("random_invertible")

#dso$optim_solution(iterations = 5000, config=optim_config(coef_hinge_H = 1, coef_hinge_W = 1,solution_balancing_threshold = 1e26,method = "positivity"))


solution <- dso$finalize_solution()
dso$plot_projected(
  "zero_distance",
  "zero_distance",
  with_solution = TRUE,
  use_dims = list(2:3)
)
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

dso$plot_negative_basis_change()
dso$plot_negative_proportions_change()






#Data preparation
n_ct <- 3
dataset <- "GSE19830"
gse <- getGEO(dataset, AnnotGPL = T)
gse <- gse[[1]]
data_raw <- Biobase::exprs(gse)
#linearize the data!
data_raw <- linearize_dataset(data_raw)

## Map probe IDs to genes
# probes which are present in data
probe_mapping <- biomaRt::select(rat2302.db::rat2302.db, rownames(data_raw), c("SYMBOL", "GENETYPE"))
probe_mapping <- probe_mapping[probe_mapping$PROBEID %in% rownames(data_raw),]
# not empty result gene id (now we  many-to-many in in genes to probes)
probe_mapping <- probe_mapping[!is.na(probe_mapping$SYMBOL),]
# This dataset is strange, 1 probe could be many genes, which is not correct (but you can try without this)
probe_mapping <- probe_mapping[!duplicated(probe_mapping$PROBEID),]
## want to collapse gene for each gene and take mean value of data
data_raw <- data_raw[probe_mapping$PROBEID, ]
rownames(data_raw) <- probe_mapping$SYMBOL
# should be unique gene name as row name after this
data_raw <- tapply(data_raw, list(row.names(data_raw)[row(data_raw)], colnames(data_raw)[col(data_raw)]), FUN = median)

# H
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

# W
components <- lapply(c(1:n_ct), function(comp_num){
  component_subset_columns <- colnames(true_proportions[,true_proportions[comp_num, ] == 1])
  component_subset <- data_raw[, component_subset_columns]
  component_vector <- as.matrix(rowMedians(component_subset))
  rownames(component_vector) <- rownames(component_subset)
  return(component_vector)
  
})
true_basis <- do.call(cbind, components)
colnames(true_basis) <- component_names

print(paste("Dim for the whole data:", toString(dim(data_raw))))
print(paste("Dim for hidden proportions (H):", toString(dim(true_proportions))))
print(paste("Dim for hidden basis (W):", toString(dim(true_basis))))

coding_gene_info = probe_mapping[probe_mapping$GENETYPE == "protein-coding",]

gene_anno <- list(
  CODING = coding_gene_info$SYMBOL
)

dso <- DualSimplexSolver$new()
dso$set_data(data_raw, gene_anno_lists = gene_anno, max_dim = 20)
dso$basic_filter(remove_true_cols_default = c(), keep_true_cols= c("CODING"))
dso$project(n_ct)

print(paste("Dim for the prefiltered data data:", toString(dim(dso$get_data()))))

dso$plot_projected("zero_distance", "zero_distance", use_dims = list(2:3))

dso$add_density_anno(genes = T)
density_values <- dso[["st"]][["data"]]@featureData@data[["density"]]
df <- data.frame(density = density_values)


# Запуск итеративной фильтрации
dso <- iterative_density_filter(dso, 
                                threshold = 1, 
                                n_cell_types = n_ct, 
                                genes = TRUE)

# Повторная проекция и визуализация
dso$project(n_ct)
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