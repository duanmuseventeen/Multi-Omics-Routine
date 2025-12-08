rm(list = ls());gc()

# Load pkgs---------------------------------------------------------------------
set.seed(1011)
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
library(harmony)
library(stringr)
# library(decontX)
library(scDblFinder)
library(DoubletFinder)
library(Augur)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(patchwork)
library(tidyr)
library(DESeq2)
library(RColorBrewer)
library(GSVA)
# library(slingshot) # Trajectory
# library(tradeSeq) # Trajectory
# library(monocle) # Trajectory
library(monocle3) # Trajectory
library("CellChat") # Communication
# library("iTALK") # Communication
# library(infercnv) # CNV
# library(copykat) # CNV
# Load data---------------------------------------------------------------------
set.seed(1011)
# options(future.globals.maxSize = 1000 * 1024^2) # 设置为1GB

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

# linux shell ------------------------------------------------------------------
for i in $(seq 32 47); do 
  mkdir GSM46795${i} 
  mv GSM46795${i}_*gz GSM46795${i}
  mv GSM46795${i}/GSM46795${i}_*barcodes.tsv.gz GSM46795${i}/barcodes.tsv.gz
  mv GSM46795${i}/GSM46795${i}_*features.tsv.gz GSM46795${i}/features.tsv.gz
  mv GSM46795${i}/GSM46795${i}_*genes.tsv.gz    GSM46795${i}/features.tsv.gz
  mv GSM46795${i}/GSM46795${i}_*matrix.mtx.gz   GSM46795${i}/matrix.mtx.gz
done
# R ----------------------------------------------------------------------------
hub        <- dir()
names(hub) <- hub
scList     <- lapply(hub, function(x){
  sce <- CreateSeuratObject(
    counts = Read10X(x),
    project = x,
    min.cells = 0,
    min.features = 0,
    assay = "RNA")
  return(sce)
})
lapply(scList, dim)
lapply(scList, function(x){ rownames(x) }) %>% Reduce(intersect, .)  -> commonsymbol
length(commonsymbol) # [1] 20283
lapply(scList, function(x){ x <- x[rownames(x) %in% commonsymbol,] })-> scList.com
lapply(scList.com, dim)

scobj      <- merge(x = scList.com[[1]], y = scList.com[-1], add.cell.ids = hub)
# raw data stat-----------------------------------------------------------------
scobj
# An object of class Seurat 
# 20283 features across 17086 samples within 1 assay 
# Active assay: RNA (20283 features, 0 variable features)
# 16 layers present: counts.GSM4679532, counts.GSM4679533, counts.GSM4679534, 
# counts.GSM4679535, counts.GSM4679536, counts.GSM4679537, counts.GSM4679538, 
# counts.GSM4679539, counts.GSM4679540, counts.GSM4679541, counts.GSM4679542, 
# counts.GSM4679543, counts.GSM4679544, counts.GSM4679545, counts.GSM4679546, 
# counts.GSM4679547

dim(scobj)
# [1] 20283 17086
# scobj <- scobj[!str_detect(rownames(scobj), "^ENSG"),]

table(scobj$orig.ident)
# GSM4679532 GSM4679533 GSM4679534 GSM4679535 GSM4679536 GSM4679537 GSM4679538 
# 585        786        837       1026        913        769       1098 
# GSM4679539 GSM4679540 GSM4679541 GSM4679542 GSM4679543 GSM4679544 GSM4679545 
# 1139        898       1570        533        745        526        272 
# GSM4679546 GSM4679547 
# 2905       2484
# QC----------------------------------------------------------------------------
setwd("2.downstream")
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

pdf("run1|QC.pdf")
VlnPlot(scobj, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(scobj, features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, features = c("percent.rp"), ncol = 1)
VlnPlot(scobj, features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
