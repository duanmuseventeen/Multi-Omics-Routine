# ST 和 scRNSseq处理中的区别：
# 1. ST中的最小单位是spot
# 2. ST中spot是多细胞混合，所以需要去卷积 或者 用配套的scRNAseq估计概率

rm(list = ls())

setwd("E:/cold2hot/3.下游处理(表达矩阵数据)/scRNAseq/")
# load pkgs --------------------------------------------------------------------
require(Seurat)
require(dplyr)
require(ggplot2)
require(stringr)
# load data --------------------------------------------------------------------
myst <- readRDS("GSE235672/GSE235672_GBM.spatial.rds")

class(myst)
# [1] "Seurat"
# attr(,"package")
# [1] "SeuratObject"

Assays(myst)
# [1] "Spatial" "SCT" 

str(myst)

View(myst)

myst@meta.data %>% dim
myst@meta.data %>% head

table(myst$SCT_snn_res.0.8, myst$seurat_clusters)
#       0    1   10   11   12    2    3    4    5    6    7    8    9
# 0  8627    0    0    0    0    0    0    0    0    0    0    0    0
# 1     0 6890    0    0    0    0    0    0    0    0    0    0    0
# 10    0    0  485    0    0    0    0    0    0    0    0    0    0
# 11    0    0    0  203    0    0    0    0    0    0    0    0    0
# 12    0    0    0    0   66    0    0    0    0    0    0    0    0
# 2     0    0    0    0    0 6005    0    0    0    0    0    0    0
# 3     0    0    0    0    0    0 5435    0    0    0    0    0    0
# 4     0    0    0    0    0    0    0 4891    0    0    0    0    0
# 5     0    0    0    0    0    0    0    0 4162    0    0    0    0
# 6     0    0    0    0    0    0    0    0    0 3290    0    0    0
# 7     0    0    0    0    0    0    0    0    0    0 2504    0    0
# 8     0    0    0    0    0    0    0    0    0    0    0 1615    0
# 9     0    0    0    0    0    0    0    0    0    0    0    0  656
# Reduction     ----------------------------------------------------------------
myst@assays$Spatial@counts %>% dim # [1] 38224 44829
myst@assays$SCT@counts %>% dim     # [1] 23700 44829

all(rownames(myst@assays$SCT@counts) %in% rownames(myst@assays$Spatial@counts)) # [1] TRUE
identical(colnames(myst@assays$SCT@counts) , colnames(myst@assays$Spatial@counts))  # [1] TRUE

myst <- NormalizeData(myst)
myst <- FindVariableFeatures(myst)
myst <- ScaleData(myst)
myst <- RunPCA(myst)
myst <- FindNeighbors(myst, dims = 1:30)
myst <- FindClusters(myst, resolution = 0.5)
myst <- RunUMAP(myst, dims = 1:30)
# Visualization ----------------------------------------------------------------
VlnPlot(myst, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(myst, features = "nCount_Spatial", ncol = 10) + 
  NoLegend() +
  theme(
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
SpatialFeaturePlot(myst, features = "HPCA") + 
  NoLegend() +
  theme(
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
SpatialFeaturePlot(myst, features = "TTR") + 
  NoLegend() +
  theme(
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  )






