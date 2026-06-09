# Noise-Robust Topological Gene Expression Deconvolution

This repository contains the implementation and experimental code for the project  
**“Development of a Noise-Robust Topological Gene Expression Deconvolution Approach”**.

The project introduces a preprocessing framework designed to improve the robustness of
**reference-free gene expression deconvolution methods**, with a particular focus on
topology-based approaches such as DualSimplex.
<img width="974" height="477" alt="image" src="https://github.com/user-attachments/assets/aca70a9a-9f95-40b5-ba65-1d01b823dc72" />


---

## 🔬 Motivation

Bulk transcriptomic deconvolution aims to recover cell-type proportions from mixed gene
expression profiles. While many existing approaches perform well under moderate noise,
their accuracy deteriorates in the presence of **extreme noise and outliers**, which are
common in real-world datasets.

This work shows that:
- Deconvolution quality is particularly sensitive to **outlier-dominated noise**, rather
  than noise variance alone.
- A topology-aware, iterative density-based filtering procedure can substantially
  improve deconvolution accuracy.
- The proposed preprocessing step is especially effective for **reference-free
  deconvolution methods**, but can also be applied prior to supervised approaches.

---

## 🧠 Method Overview

The proposed pipeline consists of three main stages:

1. **Noise and outlier modeling**  
   Synthetic benchmarks are generated using multiple noise models to simulate realistic
   distortions in gene expression space.

2. **Iterative Density Filtering**  
   An iterative procedure identifies and removes extreme outliers while preserving the
   global geometric structure of the data.

3. **Reference-Free Deconvolution**  
   The filtered expression matrix is used as input for deconvolution methods (e.g.
   DualSimplex), resulting in improved estimation of:
   - cell-type proportions  
   - basis (cell-type expression) matrices

---

## 📊 Key Results

- The filtering procedure consistently improves deconvolution accuracy across most noise
  models.
- Performance gains are especially pronounced under heavy-tailed and outlier-rich noise.
- The approach generalizes across synthetic benchmarks and real bulk RNA-seq datasets
  (e.g. omnideconv benchmarks, Hoek-purified data).

---

## 🚀 Getting Started

### Requirements

- R (≥ 4.2)
- Recommended packages:
biomaRt
GEOquery
ggplot2
ComplexHeatmap
progress
reshape
plotly
ggrastr
ggpubr
linseed
svglite
DualSimplex
dbscan
TOAST
rJava
CAM3
CellDistinguisher
NMF
hNMF
matrixStats
Biobase
dplyr
RColorBrewer
hgu133plus2.db

### Examples

```r
# Run iterative density filter on simulated data
source("scripts/simulation_iterative_density_filter.R")

# Validate on GSE19830
source("scripts/realdata_gse19830_density_filter.R")

# Validate across different reference-free methods on simulated data
source("scripts/validate_reference_free_methods_simulated.R")

# Validate on GSE64655
source("scripts/GSE64655_reference_free.R")

```

