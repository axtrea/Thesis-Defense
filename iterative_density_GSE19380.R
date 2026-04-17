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
                                     n_cell_types = 3,
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
    dso$project(n_cell_types)
    
    iter <- iter + 1
  }
  
  invisible(dso)
}

source("D:/downloads/figure_utils.R")
source("D:/downloads/setup.R")
dir_to_save_fig <- "D:/downloads/out/density_GSE19380/"
dir.create(file.path(".", dir_to_save_fig), showWarnings = F, recursive = T)

#Data preparation
n_ct <- 4
dataset <- "GSE19380"
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


# Вектор идентификаторов образцов (GSM) в правильном порядке
gsm_ids <- c(
  "GSM480943", "GSM480944", "GSM480945", "GSM480946",  # N1–N4
  "GSM480947", "GSM480948", "GSM480949", "GSM480950",  # A1–A4
  "GSM480951", "GSM480952", "GSM480953", "GSM480954",  # O1–O4
  "GSM480955", "GSM480956", "GSM480957", "GSM480958",  # M1–M4
  "GSM480959", "GSM480960", "GSM480961", "GSM480962",  # смеси NA
  "GSM480963", "GSM480964", "GSM480965", "GSM480966",  # смеси NA и NAO
  "GSM480967", "GSM480968"                              # смеси NAOM
)

# Названия клеточных типов
cell_types <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia")

# Создаём матрицу, заполненную нулями
true_proportions <- matrix(0, nrow = length(cell_types), ncol = length(gsm_ids),
                           dimnames = list(cell_types, gsm_ids))

# 1. Чистые образцы (для каждого типа ставим 1 в соответствующей строке)
# Нейроны (первые 4)
true_proportions["Neurons", gsm_ids[1:4]] <- 1
# Астроциты (следующие 4)
true_proportions["Astrocytes", gsm_ids[5:8]] <- 1
# Олигодендроциты (следующие 4)
true_proportions["Oligodendrocytes", gsm_ids[9:12]] <- 1
# Микроглия (следующие 4)
true_proportions["Microglia", gsm_ids[13:16]] <- 1

# 2. Смеси (образцы 17–26)
# 17: GSM480959 – N 0.25, A 0.75
true_proportions["Neurons", "GSM480959"] <- 0.25
true_proportions["Astrocytes", "GSM480959"] <- 0.75

# 18: GSM480960 – N 0.25, A 0.75
true_proportions["Neurons", "GSM480960"] <- 0.25
true_proportions["Astrocytes", "GSM480960"] <- 0.75

# 19: GSM480961 – N 0.5, A 0.5
true_proportions["Neurons", "GSM480961"] <- 0.5
true_proportions["Astrocytes", "GSM480961"] <- 0.5

# 20: GSM480962 – N 0.5, A 0.5
true_proportions["Neurons", "GSM480962"] <- 0.5
true_proportions["Astrocytes", "GSM480962"] <- 0.5

# 21: GSM480963 – N 0.75, A 0.25
true_proportions["Neurons", "GSM480963"] <- 0.75
true_proportions["Astrocytes", "GSM480963"] <- 0.25

# 22: GSM480964 – N 0.75, A 0.25
true_proportions["Neurons", "GSM480964"] <- 0.75
true_proportions["Astrocytes", "GSM480964"] <- 0.25

# 23: GSM480965 – N 0.5, A 0.25, O 0.25
true_proportions["Neurons", "GSM480965"] <- 0.5
true_proportions["Astrocytes", "GSM480965"] <- 0.25
true_proportions["Oligodendrocytes", "GSM480965"] <- 0.25

# 24: GSM480966 – N 0.5, A 0.25, O 0.25
true_proportions["Neurons", "GSM480966"] <- 0.5
true_proportions["Astrocytes", "GSM480966"] <- 0.25
true_proportions["Oligodendrocytes", "GSM480966"] <- 0.25

# 25: GSM480967 – N 0.5, A 0.2, O 0.2, M 0.1
true_proportions["Neurons", "GSM480967"] <- 0.5
true_proportions["Astrocytes", "GSM480967"] <- 0.2
true_proportions["Oligodendrocytes", "GSM480967"] <- 0.2
true_proportions["Microglia", "GSM480967"] <- 0.1

# 26: GSM480968 – N 0.5, A 0.2, O 0.2, M 0.1
true_proportions["Neurons", "GSM480968"] <- 0.5
true_proportions["Astrocytes", "GSM480968"] <- 0.2
true_proportions["Oligodendrocytes", "GSM480968"] <- 0.2
true_proportions["Microglia", "GSM480968"] <- 0.1

# Проверим, что сумма по каждому столбцу равна 1 (с учётом погрешностей)
colSums(true_proportions)





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
  threshold = 1,
  n_cell_types = 4,
  genes = T, 
  density_radius = NULL, 
  max_iter = 100,
  verbose = TRUE
)

dso$project(4)
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


ptp <- coerce_pred_true_props(cur_H, true_proportions)
plot_ptp_scatter(ptp)
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



