# ST 和 scRNSseq处理中的区别：
# 1. ST中的最小单位是spot
# 2. ST中spot是多细胞混合，所以需要去卷积 或者 用配套的scRNAseq估计概率

rm(list = ls())

setwd("E:/cold2hot/3.下游处理(表达矩阵数据)/scRNAseq/GSE235672/RAW/")
# load pkgs --------------------------------------------------------------------
require(Seurat)
require(dplyr)
require(ggplot2)
require(stringr)
# load data --------------------------------------------------------------------
st.list <- lapply(dir(), function(path){
  res <- Load10X_Spatial(
    data.dir = path,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = paste0("slice", path)
  )
})
names(st.list) <- dir()

mystobj <- merge(x = st.list[[1]], y = st.list[-1], add.cell.ids = dir())

class(mystobj)
# [1] "Seurat"
# attr(,"package")
# [1] "SeuratObject"

dim(mystobj)
# [1] 38224 30141
# Reduction     ----------------------------------------------------------------
myst <- mystobj

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






