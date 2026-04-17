library('GEOquery')
library(biomaRt)
library(hgu133plus2.db)
library(ComplexHeatmap)
library(progress)
library(reshape)
library(plotly)
library(ggrastr)
library(ggpubr)
library(linseed)
library(matrixStats)
library(ComplexHeatmap)
library(svglite)
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
                                     genes = FALSE, 
                                     density_radius = NULL,
                                     max_iter = 100,
                                     verbose = TRUE) {
  
  
  iter <- 1
  repeat {
    # 1. Обновить аннотацию плотности на текущих данных
    dso$add_density_anno(radius = density_radius, genes = genes)
    
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
    dso$project(3)
    
    iter <- iter + 1
  }
  
  invisible(dso)
}

source("D:/downloads/figure_utils.R")
source("D:/downloads/setup.R")
dir_to_save_fig <- "D:/downloads/out/density_GSE19830/"
dir.create(file.path(".", dir_to_save_fig), showWarnings = F, recursive = T)

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

dso$plot_svd(dims=1:15)
dso$plot_projected("zero_distance", "zero_distance", use_dims = list(2:3))

dso$add_density_anno(genes = T)
density_values <- dso[["st"]][["data"]]@featureData@data[["density"]]
df <- data.frame(density = density_values)

p <- ggplot(df, aes(x = density)) +
  geom_histogram(fill = "lightblue", color = "white", bins = 200) +
  labs(title = "Gene density distribution",
       x = "Density",
       y = "Freq.") +
  theme_minimal()

# Добавить линию порога, если нужно
threshold <- 100
p <- p + geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1)

print(p)

dummy_threshold <- 1
data <- dso$get_data()
anno <- get_anno(data)
anno$PASS_FILTER <- FALSE
anno[anno$log_mad > dummy_threshold,]$PASS_FILTER <- TRUE
data <- set_anno(anno, data)
plot_feature(data, feature = "log_mad", col_by ='PASS_FILTER')

dso$basic_filter(log_mad_gt = dummy_threshold)
dso$project(3)

dso <- iterative_density_filter(
  dso, 
  threshold = 100, 
  genes = T, 
  density_radius = NULL, 
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


markers  <- dso$get_marker_genes()

to_plot <- data_raw[
  unlist(markers),
  apply(
    true_proportions,
    2,
    function(x) any(x == 1)
  )
]

dso$plot_projected(rownames(to_plot), use_dims = (2:3))