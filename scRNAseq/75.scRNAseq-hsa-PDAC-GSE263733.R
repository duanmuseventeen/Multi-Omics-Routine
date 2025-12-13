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
setwd("GSE263733")
# options(future.globals.maxSize = 1000 * 1024^2) # 设置为1GB

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

mtx <- data.table::fread("GSE263733_Raw_counts.txt")

dim(mtx) # [1] 33694 69205
sum(duplicated(mtx$V1)) # [1] 0

mtx <- as.data.frame(mtx)
rownames(mtx) <- mtx$V1
mtx <- mtx[,-1]

mtx <- mtx %>% 
  as.matrix %>% 
  Matrix(sparse = TRUE)

table(str_remove_all(colnames(mtx), "_[ACTG]+$"))
# PB2032-PAAD-T     PB2032-PAAD-T2      PB2151-PAAD-N      PB2151-PAAD-T 
# 933               1056               2773               3710 
# PB2151-PAAD-T1 PB2155-Livermeta-T      PB2155-PAAD-N      PB2155-PAAD-T 
# 2510                460                293               1973 
# PB2191-Livermeta-T      PB2191-PAAD-N      PB2191-PAAD-T      PB2203-PAAD-N 
# 373               2567               1979                372 
# PB2203-PAAD-T      PB2218-PAAD-T      PB2219-PAAD-T      PB2256-PAAD-T 
# 651                688               1772               2423 
# PB2264-Livermeta-T      PB2264-PAAD-T      PB2265-PAAD-T      PB2266-PAAD-T 
# 3066               1121               1627               2244 
# PB2268-PAAD-T      PB2281-PAAD-T      PB2286-PAAD-T      PB2287-PAAD-N 
# 1702               2833               2344                644 
# PB2287-PAAD-T PB2311-Livermeta-T      PB2311-PAAD-T      PB2341-PAAD-T 
# 2451               2596                253               3895 
# PB2349-Livermeta-T      PB2349-PAAD-T      PB2366-PAAD-T PB2409-Livermeta-T 
# 1973               2713               2857               5380 
# PB2409-PAAD-T PB2410-Livermeta-T      PB2410-PAAD-T 
# 1953               2559               2460

metadata <- data.frame(
  orig.ident = str_remove_all(colnames(mtx), "_[ACTG]+$"),
  row.names = colnames(mtx)
) %>% mutate(sample = str_remove(orig.ident, "-.*"),
             tissue = str_extract(orig.ident, "PAAD|Livermeta"))

table(metadata$orig.ident, metadata$sample)
# PB2032 PB2151 PB2155 PB2191 PB2203 PB2218 PB2219 PB2256
# PB2032-PAAD-T         933      0      0      0      0      0      0      0
# PB2032-PAAD-T2       1056      0      0      0      0      0      0      0
# PB2151-PAAD-N           0   2773      0      0      0      0      0      0
# PB2151-PAAD-T           0   3710      0      0      0      0      0      0
# PB2151-PAAD-T1          0   2510      0      0      0      0      0      0
# PB2155-Livermeta-T      0      0    460      0      0      0      0      0
# PB2155-PAAD-N           0      0    293      0      0      0      0      0
# PB2155-PAAD-T           0      0   1973      0      0      0      0      0
# PB2191-Livermeta-T      0      0      0    373      0      0      0      0
# PB2191-PAAD-N           0      0      0   2567      0      0      0      0
# PB2191-PAAD-T           0      0      0   1979      0      0      0      0
# PB2203-PAAD-N           0      0      0      0    372      0      0      0
# PB2203-PAAD-T           0      0      0      0    651      0      0      0
# PB2218-PAAD-T           0      0      0      0      0    688      0      0
# PB2219-PAAD-T           0      0      0      0      0      0   1772      0
# PB2256-PAAD-T           0      0      0      0      0      0      0   2423
# PB2264-Livermeta-T      0      0      0      0      0      0      0      0
# PB2264-PAAD-T           0      0      0      0      0      0      0      0
# PB2265-PAAD-T           0      0      0      0      0      0      0      0
# PB2266-PAAD-T           0      0      0      0      0      0      0      0
# PB2268-PAAD-T           0      0      0      0      0      0      0      0
# PB2281-PAAD-T           0      0      0      0      0      0      0      0
# PB2286-PAAD-T           0      0      0      0      0      0      0      0
# PB2287-PAAD-N           0      0      0      0      0      0      0      0
# PB2287-PAAD-T           0      0      0      0      0      0      0      0
# PB2311-Livermeta-T      0      0      0      0      0      0      0      0
# PB2311-PAAD-T           0      0      0      0      0      0      0      0
# PB2341-PAAD-T           0      0      0      0      0      0      0      0
# PB2349-Livermeta-T      0      0      0      0      0      0      0      0
# PB2349-PAAD-T           0      0      0      0      0      0      0      0
# PB2366-PAAD-T           0      0      0      0      0      0      0      0
# PB2409-Livermeta-T      0      0      0      0      0      0      0      0
# PB2409-PAAD-T           0      0      0      0      0      0      0      0
# PB2410-Livermeta-T      0      0      0      0      0      0      0      0
# PB2410-PAAD-T           0      0      0      0      0      0      0      0
# 
# PB2264 PB2265 PB2266 PB2268 PB2281 PB2286 PB2287 PB2311
# PB2032-PAAD-T           0      0      0      0      0      0      0      0
# PB2032-PAAD-T2          0      0      0      0      0      0      0      0
# PB2151-PAAD-N           0      0      0      0      0      0      0      0
# PB2151-PAAD-T           0      0      0      0      0      0      0      0
# PB2151-PAAD-T1          0      0      0      0      0      0      0      0
# PB2155-Livermeta-T      0      0      0      0      0      0      0      0
# PB2155-PAAD-N           0      0      0      0      0      0      0      0
# PB2155-PAAD-T           0      0      0      0      0      0      0      0
# PB2191-Livermeta-T      0      0      0      0      0      0      0      0
# PB2191-PAAD-N           0      0      0      0      0      0      0      0
# PB2191-PAAD-T           0      0      0      0      0      0      0      0
# PB2203-PAAD-N           0      0      0      0      0      0      0      0
# PB2203-PAAD-T           0      0      0      0      0      0      0      0
# PB2218-PAAD-T           0      0      0      0      0      0      0      0
# PB2219-PAAD-T           0      0      0      0      0      0      0      0
# PB2256-PAAD-T           0      0      0      0      0      0      0      0
# PB2264-Livermeta-T   3066      0      0      0      0      0      0      0
# PB2264-PAAD-T        1121      0      0      0      0      0      0      0
# PB2265-PAAD-T           0   1627      0      0      0      0      0      0
# PB2266-PAAD-T           0      0   2244      0      0      0      0      0
# PB2268-PAAD-T           0      0      0   1702      0      0      0      0
# PB2281-PAAD-T           0      0      0      0   2833      0      0      0
# PB2286-PAAD-T           0      0      0      0      0   2344      0      0
# PB2287-PAAD-N           0      0      0      0      0      0    644      0
# PB2287-PAAD-T           0      0      0      0      0      0   2451      0
# PB2311-Livermeta-T      0      0      0      0      0      0      0   2596
# PB2311-PAAD-T           0      0      0      0      0      0      0    253
# PB2341-PAAD-T           0      0      0      0      0      0      0      0
# PB2349-Livermeta-T      0      0      0      0      0      0      0      0
# PB2349-PAAD-T           0      0      0      0      0      0      0      0
# PB2366-PAAD-T           0      0      0      0      0      0      0      0
# PB2409-Livermeta-T      0      0      0      0      0      0      0      0
# PB2409-PAAD-T           0      0      0      0      0      0      0      0
# PB2410-Livermeta-T      0      0      0      0      0      0      0      0
# PB2410-PAAD-T           0      0      0      0      0      0      0      0
# 
# PB2341 PB2349 PB2366 PB2409 PB2410
# PB2032-PAAD-T           0      0      0      0      0
# PB2032-PAAD-T2          0      0      0      0      0
# PB2151-PAAD-N           0      0      0      0      0
# PB2151-PAAD-T           0      0      0      0      0
# PB2151-PAAD-T1          0      0      0      0      0
# PB2155-Livermeta-T      0      0      0      0      0
# PB2155-PAAD-N           0      0      0      0      0
# PB2155-PAAD-T           0      0      0      0      0
# PB2191-Livermeta-T      0      0      0      0      0
# PB2191-PAAD-N           0      0      0      0      0
# PB2191-PAAD-T           0      0      0      0      0
# PB2203-PAAD-N           0      0      0      0      0
# PB2203-PAAD-T           0      0      0      0      0
# PB2218-PAAD-T           0      0      0      0      0
# PB2219-PAAD-T           0      0      0      0      0
# PB2256-PAAD-T           0      0      0      0      0
# PB2264-Livermeta-T      0      0      0      0      0
# PB2264-PAAD-T           0      0      0      0      0
# PB2265-PAAD-T           0      0      0      0      0
# PB2266-PAAD-T           0      0      0      0      0
# PB2268-PAAD-T           0      0      0      0      0
# PB2281-PAAD-T           0      0      0      0      0
# PB2286-PAAD-T           0      0      0      0      0
# PB2287-PAAD-N           0      0      0      0      0
# PB2287-PAAD-T           0      0      0      0      0
# PB2311-Livermeta-T      0      0      0      0      0
# PB2311-PAAD-T           0      0      0      0      0
# PB2341-PAAD-T        3895      0      0      0      0
# PB2349-Livermeta-T      0   1973      0      0      0
# PB2349-PAAD-T           0   2713      0      0      0
# PB2366-PAAD-T           0      0   2857      0      0
# PB2409-Livermeta-T      0      0      0   5380      0
# PB2409-PAAD-T           0      0      0   1953      0
# PB2410-Livermeta-T      0      0      0      0   2559
# PB2410-PAAD-T           0      0      0      0   2460
# Create Seurat object----------------------------------------------------------
scobj <- CreateSeuratObject(
  counts = mtx, 
  meta.data = metadata, 
  project = "GSE263733")

dim(scobj) # [1] 33694 69204

scobj <- scobj[!str_detect(rownames(scobj), "^A[CFJLP][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^B[X][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^ENSG"),]

dim(scobj) # [1] 31755 69204

table(scobj$orig.ident)
# PB2032-PAAD-T     PB2032-PAAD-T2      PB2151-PAAD-N      PB2151-PAAD-T 
# 933               1056               2773               3710 
# PB2151-PAAD-T1 PB2155-Livermeta-T      PB2155-PAAD-N      PB2155-PAAD-T 
# 2510                460                293               1973 
# PB2191-Livermeta-T      PB2191-PAAD-N      PB2191-PAAD-T      PB2203-PAAD-N 
# 373               2567               1979                372 
# PB2203-PAAD-T      PB2218-PAAD-T      PB2219-PAAD-T      PB2256-PAAD-T 
# 651                688               1772               2423 
# PB2264-Livermeta-T      PB2264-PAAD-T      PB2265-PAAD-T      PB2266-PAAD-T 
# 3066               1121               1627               2244 
# PB2268-PAAD-T      PB2281-PAAD-T      PB2286-PAAD-T      PB2287-PAAD-N 
# 1702               2833               2344                644 
# PB2287-PAAD-T PB2311-Livermeta-T      PB2311-PAAD-T      PB2341-PAAD-T 
# 2451               2596                253               3895 
# PB2349-Livermeta-T      PB2349-PAAD-T      PB2366-PAAD-T PB2409-Livermeta-T 
# 1973               2713               2857               5380 
# PB2409-PAAD-T PB2410-Livermeta-T      PB2410-PAAD-T 
# 1953               2559               2460 
# QC----------------------------------------------------------------------------
setwd("GSE263733/2.downstream")

scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

pdf("run1|QC.pdf")
VlnPlot(scobj, features = c("nFeature_RNA"), ncol = 1) + NoLegend()
VlnPlot(scobj, features = c("nCount_RNA"), ncol = 1) + NoLegend()
VlnPlot(scobj, features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20)) + NoLegend()
VlnPlot(scobj, features = c("percent.rp"), ncol = 1) + NoLegend()
VlnPlot(scobj, features = c("percent.hb"), ncol = 1) + NoLegend()
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 500 &
    nFeature_RNA < 7500 & 
    nCount_RNA > 500 &
    nCount_RNA < 30000 &
    percent.mt < 20 &
    percent.rp < 40 &
    percent.hb < 1)  

pdf("run1|QC R2.pdf")
VlnPlot(scobj, features = c("nFeature_RNA"), ncol = 1) + NoLegend()
VlnPlot(scobj, features = c("nCount_RNA"), ncol = 1) + NoLegend()
VlnPlot(scobj, features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20)) + NoLegend()
VlnPlot(scobj, features = c("percent.rp"), ncol = 1) + NoLegend()
VlnPlot(scobj, features = c("percent.hb"), ncol = 1) + NoLegend()
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

table(scobj$orig.ident)
# PB2032-PAAD-T     PB2032-PAAD-T2      PB2151-PAAD-N      PB2151-PAAD-T 
# 818                925               2362               3152 
# PB2151-PAAD-T1 PB2155-Livermeta-T      PB2155-PAAD-N      PB2155-PAAD-T 
# 2061                298                213               1606 
# PB2191-Livermeta-T      PB2191-PAAD-N      PB2191-PAAD-T      PB2203-PAAD-N 
# 220               2312               1660                340 
# PB2203-PAAD-T      PB2218-PAAD-T      PB2219-PAAD-T      PB2256-PAAD-T 
# 577                529               1565               2106 
# PB2264-Livermeta-T      PB2264-PAAD-T      PB2265-PAAD-T      PB2266-PAAD-T 
# 2344                836               1122               1831 
# PB2268-PAAD-T      PB2281-PAAD-T      PB2286-PAAD-T      PB2287-PAAD-N 
# 1281               2426               1830                536 
# PB2287-PAAD-T PB2311-Livermeta-T      PB2311-PAAD-T      PB2341-PAAD-T 
# 2110               2137                 78               2884 
# PB2349-Livermeta-T      PB2349-PAAD-T      PB2366-PAAD-T PB2409-Livermeta-T 
# 1587               2496               2385               4858 
# PB2409-PAAD-T PB2410-Livermeta-T      PB2410-PAAD-T 
# 1506               1998               2101 

dim(scobj) # [1] 31755 57090
# blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(scobj)) # [1] 240
scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj)
# [1] 31516 57090

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "run1|scobj.Rdata")
# Nomalization & harmony & cluster----------------------------------------------
scobj.harmony <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "sample"),
    reduction.use = "pca",
    reduction.save = "harmony") %>% 
  JoinLayers(assay = "RNA")

scobj.harmony.20 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.25 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.30 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

pdf("run1|nf2000_h_clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony, reduction = "pca", ndims = 50)
clustree(scobj.harmony.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.30, prefix = "RNA_snn_res.")
dev.off()

scobj.harmony <- scobj.harmony.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Visualization-----------------------------------------------------------------
pdf("run1|nf2000_h_pc30 umap.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.1", reduction = "umap", label = T)
FeaturePlot(scobj.harmony, features = "percent.mt")
dev.off()
# Cell Cycle--------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.harmony <- CellCycleScoring(scobj.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scobj.harmony <- RunPCA(scobj.harmony, features = c(s.genes, g2m.genes), reduction.name = "cellcycle")

pdf("run1|nf2000_h_pc30 umap cellcycle.pdf")
RidgePlot(scobj.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
DimPlot(scobj.harmony, reduction = "cellcycle")
DimPlot(scobj.harmony, reduction = "umap")
dev.off()
# doublets----------------------------------------------------------------------
save(scobj.harmony, file = "run1|nf2000_h_pc30 scobj.harmony.Rdata")

scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list

for (i in 1:length(scobj.harmony.split)) {
  # pK Identification (ground-truth) -------------------------------------------
  sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:30)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
  ## Homotypic Doublet Proportion Estimate -------------------------------------
  homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]@meta.data$RNA_snn_res.0.1)  
  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi <- round(0.05 * nrow(scobj.harmony.split[[i]]@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies ----------------
  # scobj.harmony.split[[i]] <- doubletFinder(scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/228
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/225#issuecomment-2786505997
  scobj.harmony.split[[i]] <- doubletFinder(
    scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, 
    nExp = nExp_poi.adj, sct = FALSE)
}

save(scobj.harmony.split, file = "run1|nf2000_h_pc30 scobj.harmony.split.Rdata")


Singlet <- c()
for (i in 1:length(scobj.harmony.split)) {
  Singlet <- c(Singlet, 
               rownames(scobj.harmony.split[[i]]@meta.data) [scobj.harmony.split[[i]]@meta.data$DF.classifications_0.25 == "Singlet"])
}
for (i in 1:length(scobj.harmony.split)) {
  all(scobj.harmony.split[[i]]@meta.data[rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet,   25] == "Singlet") %>% print
  all(scobj.harmony.split[[i]]@meta.data[!(rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet),25] != "Singlet") %>% print
}

scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)

# scobj.h.dblfinder <- subset(scobj.harmony, subset = id %in% Singlet)
scobj.harmony@meta.data$dblfinder <- NA
scobj.harmony@meta.data$dblfinder <- "doublet"
scobj.harmony@meta.data$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"

dim(scobj.harmony) # [1] 31516 57090
table(scobj.harmony@meta.data$dblfinder)
# doublet singlet 
# 1878   55212
# Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1|scobj4run2.Rdata")
## regcc========================================================================
# 1000 regccmt----
set.seed(1011)
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
library(harmony)
library(stringr)

load("GSE263733/2.downstream/run1|scobj4run2.Rdata")
scobj.harmony.1000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "sample"),
    reduction.use = "pca",
    reduction.save = "harmony") %>% 
  JoinLayers(assay = "RNA")
scobj.harmony.1000.20 <- scobj.harmony.1000 %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.1000.25 <- scobj.harmony.1000 %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.1000.30 <- scobj.harmony.1000 %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
pdf("run2|nf1000_h_regccmt_clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony.1000, 
     scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
     file = "run2|nf1000_h_regccmt scobj.Rdata")
print("--->Finished<---")

# 1500 regccmt----
set.seed(1011)
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
library(harmony)
library(stringr)

load("GSE263733/2.downstream/run1|scobj4run2.Rdata")
scobj.harmony.1500 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "sample"),
    reduction.use = "pca",
    reduction.save = "harmony") %>% 
  JoinLayers(assay = "RNA")

scobj.harmony.1500.20 <- scobj.harmony.1500 %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.1500.25 <- scobj.harmony.1500 %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.1500.30 <- scobj.harmony.1500 %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
pdf("run2|nf1500_h_regccmt_clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony.1500, 
     scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
     file = "run2|nf1500_h_regccmt scobj.Rdata")
print("--->Finished<---")

# 2000 regccmt----
set.seed(1011)
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
library(harmony)
library(stringr)

load("GSE263733/2.downstream/run1|scobj4run2.Rdata")
scobj.harmony.2000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "sample"),
    reduction.use = "pca",
    reduction.save = "harmony") %>% 
  JoinLayers(assay = "RNA")
scobj.harmony.2000.20 <- scobj.harmony.2000 %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.2000.25 <- scobj.harmony.2000 %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.2000.30 <- scobj.harmony.2000 %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
pdf("run2|nf2000_h_regccmt_clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony.2000, 
     scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
     file = "run2|nf2000_h_regccmt scobj.Rdata")
print("--->Finished<---")

# Run Umap----
# with regression for cc and mt
# 1000 30 0.3

scobj.harmony <- scobj.harmony.1000.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Annotation====================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.2")

islets <- c("INS1", "PAX6")
acinar <- c("CELA1", "CPA1")
ductal <- c("DCDC2A", "PROX1","HNF1B","CFTR")
adipocyte   <- c("PLIN1")
meta_neoplastic <- c("ONECUT2")
schwann<- c("SOX10", "S100B")
endocrine   <- c("CHGA","GCG","IAPP","SST","PPY")
neural <- c("TAC1", "CHAT")
epithelial  <- c("EPCAM", "SFN", "KRT5", "KRT14") #, "SPRR3"
endothelial <- c("VWF", "PECAM1", "ENG", "CDH5","FLT1","FLT4") #, "CCL14"
mesothelial <- c("PDPN", "MSLN", "UPK3B")
fibroblast  <- c("FN1", "DCN", "COL1A1", "COL1A2") #, "COL3A1", "COL6A1"
SMC    <- c("TAGLN", "CNN1", "PRKG1", "FOXP2")
pericyte <- c("RGS5", "MCAM", "ACTA2", "MYH11")
T_cell <- c("CD2", "CD3D", "CD3E", "CD3G") # 0 4 10 14 
B_cell <-	c("CD19", "CD79A", "MS4A1", "CD79B") # 12
plasma <-	c("JCHAIN", "MZB1", "IGHG1", "SDC1") # 8 16 , "CD79A"
# monocyte and macrophage
myeloid<- c("CD68", "CD163", "LYZ", "CD14","FCGR3A","C1QA","C1QB","CST3") # , "IL3RA", "LAMP3", "CLEC4C", "GCA"
Neutrophil <- c("CSF3R","CXCL8","G0S2","IFITM2") # "FCGR3B","FPR1","BASP1","CXCR1","CXCR2","S100A11"
DC     <- c('CCR7',	'CLEC9A', 'CD1C',	'IRF7',	'LILRA4')
Mast   <-	c('TPSAB1','CPA3','HPGDS','KIT')# 'VWA5A','SLC18A2','HDC','CAPG','RGS13','IL1RL1','FOSB','GATA2'

p_dotplot <- DotPlot(scobj.harmony, features = c(
  islets,acinar,ductal,adipocyte,meta_neoplastic,schwann,endocrine,neural,
  epithelial, endothelial, mesothelial,fibroblast, SMC, pericyte,
  T_cell, B_cell,plasma, myeloid, Neutrophil, Mast #, DC
) %>% unique, 
group.by = "RNA_snn_res.0.3") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p_dotplot, filename = "DotPlot4Anno.pdf", width=12, height=15, units="in")

C7 <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.3", ident.1 = "7")
C9 <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.3", ident.1 = "9")
C19 <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.3", ident.1 = "19")

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(0)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(1)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(2)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(3)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(4)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(5)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(6)] <- "B cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(7)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(8)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(9)] <- "Doublets"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(10)]<- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(11)]<- "Epithelial" # "Dcutal"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(12)]<- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(13)]<- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(14)]<- "B cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(15)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(16)]<- "Mast"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(17)]<- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(18)]<- "Acinar"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(19)]<- "Doublets"

# save(scobj.harmony, file = "scobj.harmony(singlet_nf1000_h_regccmt_pc30_res0.3)anno.Rdata")
load("scobj.harmony(singlet_nf1000_h_regccmt_pc25_res0.2)anno.Rdata")

table(scobj.harmony$RNA_snn_res.0.3)
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 11587  9212  9211  6173  5117  4484  2215  2003  1797  1093   891   686   653 
# 13    14    15    16    17    18    19 
# 556   370   243   242   215   187   155

table(scobj.harmony$cell_type)
# Acinar      B cell    Doublets Endothelial  Epithelial        Mast 
# 187        2585        1248         556       12143         242 
# Mesenchyme     Myeloid      T cell 
# 891       10254       28984 
