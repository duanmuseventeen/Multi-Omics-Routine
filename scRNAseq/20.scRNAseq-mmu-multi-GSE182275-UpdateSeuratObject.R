# Tissue-resident memory CD8+ T cells possess unique transcriptional, epigenetic and functional adaptations to different tissue environments
# PMID: 35761084

# Reference
# https://blog.csdn.net/weixin_56845253/article/details/130367503

# Load pkgs---------------------------------------------------------------------
library(Seurat)
library(dplyr)
# Load data---------------------------------------------------------------------
GSE182275 <- readRDS("GSE182275_seurat_object.rds") %>% 
  UpdateSeuratObject() 
# Visualization-----------------------------------------------------------------
DimPlot(GSE182275, group.by = "tissue")
DimPlot(GSE182275, group.by = "sample")
DimPlot(GSE182275, group.by = "batch")
