setwd("GSE141445/")

rm(list = ls())
set.seed(1011)
# load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
# library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
library(decontX)
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
require(stringr)
require(Matrix)
# load data---------------------------------------------------------------------
mtx <- data.table::fread("GSM4203181_data.matrix.txt.gz")
# raw.mtx <- data.table::fread("GSM4203181_data.raw.matrix.txt.gz")

dim(mtx)
# [1] 25044 36425
# dim(raw.mtx)
# # [1] 25044 36425

mtx <- as.data.frame(mtx)
rownames(mtx) <- mtx$V1
mtx <- mtx[,-1]

mtx <- mtx %>% 
  as.matrix %>% 
  Matrix(sparse = TRUE)

table(str_remove_all(colnames(mtx), "^[ACTG]+-"))
#   1   10   11   12   13    2    3    4    5    6    7    8    9 
# 1554  679 2824 5780 1277 5302 4982 4077 1339 1586 1534 3032 2458 

metadata <- data.frame(
  orig.ident = paste0("P", str_remove_all(colnames(mtx), "^[ACTG]+-")),
  row.names = colnames(mtx))
# Create Seurat object----------------------------------------------------------
scobj <- CreateSeuratObject(
  counts = mtx, 
  meta.data = metadata, 
  project = "GSE141445")

scobj
# An object of class Seurat 
# 25044 features across 36424 samples within 1 assay 
# Active assay: RNA (25044 features, 0 variable features)
# 1 layer present: counts

sum(str_detect(rownames(scobj), "^ENSG"))
# [1] 0

table(scobj$orig.ident)
# P1  P10  P11  P12  P13   P2   P3   P4   P5   P6   P7   P8   P9 
# 1554  679 2824 5780 1277 5302 4982 4077 1339 1586 1534 3032 2458 
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

pdf("run1=QC.pdf")
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA"), ncol = 1)
VlnPlot(scobj, group.by = "orig.ident", features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.rp"), ncol = 1) + scale_y_continuous(breaks = c(10,20,30,40))
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 200 &
    # nFeature_RNA < 5000 & 
    nCount_RNA > 200 &
    nCount_RNA < 12000 &
    percent.mt < 10 &
    percent.rp < 20 &
    percent.hb < 0.5)  

pdf("run1=QC2.pdf")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(scobj, features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, features = c("percent.rp"), ncol = 1)
VlnPlot(scobj, features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

table(scobj$orig.ident)
# P1  P10  P11  P12  P13   P2   P3   P4   P5   P6   P7   P8   P9 
# 1405  554 2765 5780 1244 5095 4908 3473 1139 1029 1360 3030 1937 

# dim(scobj)
# [1] 25044 33719
# blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("C:/D/R project/Multi-Omics-Routine/scRNAseq/blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(scobj))
# [1] 222

scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj)
# [1] 24823 33719

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "run1=scobj(GSE141445).Rdata")
# Nomalization & harmony & cluster----------------------------------------------
scobj.harmony <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = "orig.ident",
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

pdf("vst.2000.h clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony, reduction = "pca", ndims = 50)
clustree(scobj.harmony.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.30, prefix = "RNA_snn_res.")
dev.off()

scobj.harmony <- scobj.harmony.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Visualization-----------------------------------------------------------------
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.2")
pdf("run1=nf2000_h clustree 20 25 30.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.2", reduction = "umap", label = T)
dev.off()
# Cell Cycle--------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.harmony <- CellCycleScoring(scobj.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf("run1=cell cycle.pdf")
RidgePlot(scobj.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
scobj.harmony <- RunPCA(scobj.harmony, features = c(s.genes, g2m.genes), reduction.name = "cellcycle")
DimPlot(scobj.harmony, reduction = "cellcycle")
DimPlot(scobj.harmony, reduction = "umap")
dev.off()
# doublets----------------------------------------------------------------------
save(scobj.harmony, file = "run1=scobj.harmony(GSE141445).Rdata")

if(TRUE){
  require(DoubletFinder)
  scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list
  
  for (i in 1:length(scobj.harmony.split)) {
    # pK Identification (ground-truth) -------------------------------------------
    sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:30)
    sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
    ## Homotypic Doublet Proportion Estimate -------------------------------------
    homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]@meta.data$RNA_snn_res.0.2)  
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
  
  save(scobj.harmony.split, file = "run1=scobj.harmony.split(GSE141445).Rdata")
  
  Singlet <- c(
    rownames(scobj.harmony.split[[1]]@meta.data) [scobj.harmony.split[[1]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[2]]@meta.data) [scobj.harmony.split[[2]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[3]]@meta.data) [scobj.harmony.split[[3]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[4]]@meta.data) [scobj.harmony.split[[4]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[5]]@meta.data) [scobj.harmony.split[[5]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[6]]@meta.data) [scobj.harmony.split[[6]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[7]]@meta.data) [scobj.harmony.split[[7]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[8]]@meta.data) [scobj.harmony.split[[8]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[9]]@meta.data) [scobj.harmony.split[[9]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[10]]@meta.data)[scobj.harmony.split[[10]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[11]]@meta.data)[scobj.harmony.split[[11]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[12]]@meta.data)[scobj.harmony.split[[12]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[13]]@meta.data)[scobj.harmony.split[[13]]@meta.data$DF.classifications_0.25 == "Singlet"])
  
  scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)
  
  # scobj.h.dblfinder <- subset(scobj.harmony, subset = id %in% Singlet)
  scobj.harmony@meta.data$dblfinder <- NA
  scobj.harmony@meta.data$dblfinder <- "doublet"
  scobj.harmony@meta.data$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"
  
  dim(scobj.harmony)
  # [1] 25044 33719
  table(scobj.harmony@meta.data$dblfinder)
  # doublet singlet 
  # 1048   32671 
}
# Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1=scobj4run2(GSE141445).Rdata")

# 1000 without regress cell cycle
if(TRUE){
  scobj.harmony.1000 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
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
  
  pdf("vst.1000.h_clustree 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.1000,
       scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
       file = "scobj.harmony.1000(without regression).Rdata")
}
# 1500 without regress cell cycle
if(TRUE){
  scobj.harmony.1500 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
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
  
  pdf("vst.1500.h_clustree 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.1500,
       scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
       file = "scobj.harmony.1500(without regression).Rdata")
}
# 2000 without regress cell cycle
if(TRUE){
  scobj.harmony.2000 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
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
  
  pdf("vst.2000.h_clustree 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.2000,
       scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
       file = "scobj.harmony.2000(without regression).Rdata")
}

# 1000
if(FALSE){
  scobj.harmony.1000 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), 
              features = rownames(scobj)) %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
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
  
  pdf("vst.1000.h_clustree 20 25 30(singlet_regcc).pdf")
  ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.1000, 
       scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
       file = "scobj.harmony.1000(reg.cc).Rdata")
}
# 1500
if(FALSE){
  scobj.harmony.1500 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), 
              features = rownames(scobj)) %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
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
  
  pdf("vst.1500.h_clustree 20 25 30(singlet_regcc).pdf")
  ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.1500, 
       scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
       file = "scobj.harmony.1500(reg.cc).Rdata")
}
# 2000
if(FALSE){
  scobj.harmony.2000 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), 
              features = rownames(scobj)) %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
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
  
  pdf("vst.2000.h_clustree 20 25 30(singlet_regcc).pdf")
  ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.2000, 
       scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
       file = "scobj.harmony.2000(reg.cc).Rdata")
}

if(regression cell cycle){
  scobj.harmony <- scobj.harmony.2000.30 %>% 
    RunUMAP(reduction = "harmony", dims = 1:30)
  # Visualization=================================================================
  scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.3")
  # DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
  # DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T)
  # Annotation====================================================================
  scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.3")
  
  epithelial <- c("EPCAM", "SFN", "KRT5", "KRT14","KLK3") #, "SPRR3"
  endothelial <- c("VWF", "PECAM1", "ENG", "CDH5") #, "CCL14"
  fibroblast <- c("FN1", "DCN", "COL1A1", "COL1A2") #, "COL3A1", "COL6A1"
  SMC <- c("TAGLN", "CNN1", "PRKG1", "FOXP2")
  pericyte <- c("RGS5", "MCAM", "ACTA2", "MYH11")
  T_cell <- c("CD2", "CD3D", "CD3E", "CD3G") # 0 4 10 14 
  B_cell <-	c("CD19", "CD79A", "MS4A1", "CD79B") # 12
  plasma <-	c("JCHAIN", "MZB1", "IGHG1", "SDC1") # 8 16 , "CD79A"
  # monocyte and macrophage
  myeloid <- c("CD68", "CD163", "LYZ", "CD14","FCGR3A","C1QA","C1QB","CST3") # , "IL3RA", "LAMP3", "CLEC4C", "GCA"
  Neutrophil <- c("CSF3R","CXCL8","G0S2","IFITM2") # "FCGR3B","FPR1","BASP1","CXCR1","CXCR2","S100A11"
  DC <- c('CCR7',	'CLEC9A', 'CD1C',	'IRF7',	'LILRA4')
  Mast <-	c('TPSAB1','CPA3','HPGDS','KIT')# 'VWA5A','SLC18A2','HDC','CAPG','RGS13','IL1RL1','FOSB','GATA2'
  Neural <- c("PLP1","NRNX1","NRNX2","NRNX3")
  Proliferation <- c("TOP2A","BIRC5","MKI67","DLGAP5")
  
  pdf("Marker dotplot.pdf")
  DotPlot(scobj.harmony, features = c(
    epithelial, endothelial, fibroblast, SMC, pericyte,
    T_cell, B_cell,plasma, myeloid, Neutrophil, Mast, DC, Neural,
    Proliferation
  ) %>% unique, 
  group.by = "RNA_snn_res.0.3") + 
    scale_color_viridis() +
    labs(x = "", y = "") +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
  
  scobj.harmony$cell_type <- "Unknown"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(0)] <- "Epithelial"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(1)] <- "T cell"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(10)]<- "Proliferation"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(11)]<- "B cell"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(12)]<- "Endothelial"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(13)]<- "Mesenchyme" # Endothelial
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(14)]<- "Mesenchyme" # Fibroblast
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(15)]<- "Mesenchyme" # Fibroblast
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(2)] <- "Endothelial"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(3)] <- "Epithelial" 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(4)] <- "Mesenchyme" # Pericyte
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(5)] <- "Epithelial"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(6)] <- "Myeloid"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(7)] <- "Mast"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(8)] <- "Epithelial"
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.3 %in% c(9)] <- "Epithelial"
  
  # save(scobj.harmony, file = "scobj.harmony(singlet-reg.cc-2000-30-res_0.3) anno.Rdata")
  load("scobj.harmony(singlet-reg.cc-2000-30-res_0.3) anno.Rdata")
  
  table(scobj.harmony$RNA_snn_res.0.3)
  # 0     1    10    11    12    13    14    15     2     3     4     5     6
  # 15077  3122   651   478   414   256   154   107  2535  2281  1542  1432  1429
  # 7     8     9
  # 1410   895   888
  table(scobj.harmony$cell_type)
  # B cell   Endothelial    Epithelial          Mast    Mesenchyme       Myeloid 
  # 478          2949         20573          1410          2059          1429 
  # Proliferation        T cell 
  # 651          3122
  # Figure2=======================================================================
  require(patchwork)
  p1 <- DimPlot(scobj.harmony, reduction = "umap", 
                group.by = "cell_type", 
                label = TRUE)
  p2 <- DimPlot(scobj.harmony, reduction = "umap", 
                group.by = "RNA_snn_res.0.3", 
                label = TRUE)
  p_out1 <- p1 + p2
  ggsave(p_out1, filename = "UMAP(GSE221603).pdf", width=10, height=5, units="in")
  
  p_fea <- lapply(c(gene_symbols), 
                  function(x){FeaturePlot(scobj.harmony, features = x) + 
                      scale_colour_gradientn(
                        colours = colorRampPalette(
                          c('gray90','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410',
                            '#E52C22','#DD1D23','#C20030','#930039','#8C003A',
                            '#6F003D','#56033F'))(1000))}) %>% 
    patchwork::wrap_plots(ncol = 2, nrow = 2)
  p_vln <- VlnPlot(scobj.harmony, features = c(gene_symbols), 
                   group.by = "cell_type")
  
  ggsave(p_fea, filename = "FeaturePlot(GSE221603).pdf", width=8, height=8, units="in")
  ggsave(p_vln, filename = "VlnPlot(GSE221603).pdf", width=8, height=8, units="in")
  
  Epi_markers <- FindMarkers(scobj.harmony, group.by = "cell_type", ident.1 = "Epithelial")
  head(Epi_markers[Epi_markers$avg_log2FC > 0 & Epi_markers$p_val_adj < 0.05,], 100)
  Epi_markers[c(gene_symbols),]
}




