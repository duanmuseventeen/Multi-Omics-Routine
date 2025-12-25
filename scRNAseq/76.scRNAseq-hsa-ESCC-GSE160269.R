setwd("~/1.rawdat")
rm(list = ls());gc()

set.seed(1011)
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
# library(ROGUE) # use for determine the optimal resolution
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
# Load data---------------------------------------------------------------------
colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

meta_cd45n<- data.table::fread("GSE160269_CD45neg_cells.txt.gz")
mtx_cd45n <- data.table::fread("GSE160269_CD45neg_UMIs.txt.gz")
meta_cd45p<- data.table::fread("GSE160269_CD45pos_cells.txt.gz")
mtx_cd45p <- data.table::fread("GSE160269_CD45pos_UMIs.txt.gz")

all(meta_cd45n$ID == colnames(mtx_cd45n)[-1]) # [1] TRUE
all(meta_cd45p$ID == colnames(mtx_cd45p)[-1]) # [1] TRUE

dim(mtx_cd45n) # [1] 17012 97632
dim(mtx_cd45p) # [1]  15175 111029
# 整合基因----------------------------------------------------------------------
meta <- bind_rows(
  meta_cd45n %>% mutate(FACS = "CD45_N"),
  meta_cd45p %>% mutate(FACS = "CD45_P")
)

length(unique(c(mtx_cd45n$V1, mtx_cd45p$V1))) # [1] 17986
length(intersect(mtx_cd45n$V1, mtx_cd45p$V1)) # [1] 14201

gene.intersect <- intersect(mtx_cd45n$V1, mtx_cd45p$V1)

mtx <- mtx_cd45n %>% filter(V1 %in% gene.intersect) %>% 
  full_join(mtx_cd45p %>% filter(V1 %in% gene.intersect), by = "V1")

dim(mtx) # [1] 14201 208660

mtx <- as.data.frame(mtx)
rownames(mtx) <- mtx$V1
mtx <- mtx[,-1]
# Create Seurat Object----------------------------------------------------------
rm(mtx_cd45n)
rm(mtx_cd45p)
gc()

create.seu <- function(rawdat, meta, name){
  suppressMessages(require(Matrix))
  rawdat <- rawdat %>% 
    as.matrix %>% 
    Matrix(sparse = TRUE)
  seu.obj <- CreateSeuratObject(counts = rawdat, meta.data = meta, project = name)
  return(seu.obj)
}
scobj <- create.seu(rawdat = mtx, meta = meta, name = "GSE160269")

scobj
# An object of class Seurat
# 14201 features across 208659 samples within 1 assay
# Active assay: RNA (14201 features, 0 variable features)
# 1 layer present: counts
save(scobj, file = "run1|GSE160269.Rdata")
# QC-----------------------------------------------------------------------------
setwd("~/2.downstream/")
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

pdf("run1|QC.pdf", width = 16)
# VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), group.by = "FACS", ncol = 5)
VlnPlot(scobj, features = c("nFeature_RNA"), group.by = "sample", ncol = 1) + scale_y_continuous(breaks = c(500,1000,2000,3000,4000,5000,7000)) + NoLegend() 
VlnPlot(scobj, features = c("nCount_RNA"), group.by = "sample", ncol = 1) + scale_y_continuous(breaks = c(500,1000,10000,15000,20000,30000)) +NoLegend() 
VlnPlot(scobj, features = c("percent.mt"), group.by = "sample", ncol = 1) + scale_y_continuous(breaks = c(10,20)) + NoLegend() 
VlnPlot(scobj, features = c("percent.rp"), group.by = "sample", ncol = 1) + NoLegend() 
VlnPlot(scobj, features = c("percent.hb"), group.by = "sample", ncol = 1) + NoLegend() 
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

dim(scobj)
# [1]  14201 208659

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 1000 &
    nFeature_RNA < 5000 & 
    nCount_RNA > 1000 &
    nCount_RNA < 15000 &
    percent.mt < 10 &
    percent.rp < 40 &
    percent.hb < 1)  

dim(scobj)
# [1]  14201 157751

pdf("run1|QC R1.pdf", width = 16)
VlnPlot(scobj, features = c("nFeature_RNA"), group.by = "sample", ncol = 1) + scale_y_continuous(breaks = c(500,1000,2500,5000)) + NoLegend() 
VlnPlot(scobj, features = c("nCount_RNA"), group.by = "sample", ncol = 1) + scale_y_continuous(breaks = c(500,1000,10000,15000,20000)) + NoLegend() 
VlnPlot(scobj, features = c("percent.mt"), group.by = "sample", ncol = 1) + scale_y_continuous(breaks = c(10,20)) + NoLegend() 
VlnPlot(scobj, features = c("percent.rp"), group.by = "sample", ncol = 1) + NoLegend() 
VlnPlot(scobj, features = c("percent.hb"), group.by = "sample", ncol = 1) + NoLegend() 
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

scobj$orig.ident <- scobj$sample

table(scobj$orig.ident)
# P104T P107T  P10T  P11T P126N P126T P127N P127T P128N P128T  P12T P130N P130T 
# 2069  2108  2911  2539  4562  1247  2068  2280  3157  2617  1855  4512  1272 
# P15T  P16T  P17T  P19T   P1T  P20T  P21T  P22T  P23T  P24T  P26T  P27T  P28T 
# 2772  2374  1540  1373  1201  2029  2490  5157  2420  6034  2234  3052  2087 
# P2T  P30T  P31T  P32T  P36T  P37T  P38T  P39T  P40T  P42T  P44T  P47T  P48T 
# 2855  2295  2777  1837  2141  1216  2708  3500  4200  2926  3416  3901  2739 
# P49T   P4T  P52T  P54T  P56T  P57T   P5T  P61T  P62T  P63T  P65T  P74T  P75T 
# 1029  3164  2098  2277   577  2041  3768  2313  1710  1441  3291  1399  1875 
# P76T  P79T  P80T  P82T  P83T  P84T  P87T  P89T   P8T  P91T  P94T   P9T 
# 2312  2400  1807   816  1870  1461  2561   955  4617  1452   821  5225 

range(.Last.value)
# [1]  577 6034
# blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(scobj)) # [1] 203
scobj     <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj) # [1]  13999 157751
# Nomalization & harmony & cluster----------------------------------------------
scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)
scobj$sample   <- str_remove(scobj$sample, "[TN]$")
save(scobj, file = "run1|scobj.Rdata")

scobj.harmony <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("FACS", "orig.ident", "sample"),
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
pdf("run1|nf2000_h clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony, reduction = "pca", ndims = 50)
clustree(scobj.harmony.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony, scobj.harmony.20, scobj.harmony.25, scobj.harmony.30,
     file = "run1|nf2000_h scobj.harmony.Rdata")
print("--->Finished!<---")
# run umap----
scobj.harmony <- scobj.harmony.25 %>% 
  RunUMAP(reduction = "harmony", dims = 1:25)
# Visualization-----------------------------------------------------------------
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.2")
pdf("run1|nf2000_h_pc25 umap.pdf")
DimPlot(scobj.harmony, group.by = "FACS", reduction = "umap", label = T) + NoLegend() 
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T) + NoLegend() 
DimPlot(scobj.harmony, group.by = "sample", reduction = "umap", label = T) + NoLegend() 
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.2", reduction = "umap", label = T)
FeaturePlot(scobj.harmony, features = "percent.mt")
dev.off()
# Cell Cycle--------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.harmony <- CellCycleScoring(scobj.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf("run1|cellcycle.pdf")
RidgePlot(scobj.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
scobj.harmony <- RunPCA(scobj.harmony, features = c(s.genes, g2m.genes), reduction.name = "cellcycle")
DimPlot(scobj.harmony, reduction = "cellcycle")
DimPlot(scobj.harmony, reduction = "umap")
dev.off()
# doublets----------------------------------------------------------------------
save(scobj.harmony, file = "run1|nf2000_h_pc25 scobj.harmony.Rdata")

scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list

for (i in 1:length(scobj.harmony.split)) {
  # pK Identification (ground-truth) -------------------------------------------
  sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:25)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
  ## Homotypic Doublet Proportion Estimate -------------------------------------
  homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]@meta.data$RNA_snn_res.0.2)  
  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi <- round(0.075 * nrow(scobj.harmony.split[[i]]@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies ----------------
  # scobj.harmony.split[[i]] <- doubletFinder(scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/228
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/225#issuecomment-2786505997
  scobj.harmony.split[[i]] <- doubletFinder(
    scobj.harmony.split[[i]], PCs = 1:25, pN = 0.25, pK = pK, 
    nExp = nExp_poi.adj, sct = FALSE)
}

save(scobj.harmony.split, file = "run1|scobj.harmony.mkdbl.Rdata")

Singlet <- c()
for (i in 1:length(scobj.harmony.split)) {
  Singlet <- c(Singlet, 
               rownames(scobj.harmony.split[[i]]@meta.data) [scobj.harmony.split[[i]]@meta.data$DF.classifications_0.25 == "Singlet"])
}
for (i in 1:length(scobj.harmony.split)) {
  tmp <- scobj.harmony.split[[i]]
  all(tmp@meta.data[rownames(tmp@meta.data) %in% Singlet,   ncol(tmp@meta.data)] == "Singlet") %>% print
  all(tmp@meta.data[!(rownames(tmp@meta.data) %in% Singlet),ncol(tmp@meta.data)] != "Singlet") %>% print
}

scobj.harmony$id <- rownames(scobj.harmony@meta.data)
scobj.harmony$dblfinder <- NA
scobj.harmony$dblfinder <- "doublet"
scobj.harmony$dblfinder[scobj.harmony$id %in% Singlet] <- "singlet"

dim(scobj.harmony)
# [1]  13999 157751
table(scobj.harmony$dblfinder)
# doublet singlet 
# 9180  148571
# Re-run========================================================================
scobj[['S.Score']]   <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']]     <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1|scobj4run2.Rdata")
print("--->Finished!<---")

# without regress cell cycle====================================================
# 1000 without regress cell cycle----
scobj.harmony.1000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("FACS", "orig.ident", "sample"),
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
save(scobj.harmony.1000,
     scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
     file = "run2|nf1000_h_noreg scobj.harmony.Rdata")
  
# 1500 without regress cell cycle----
scobj.harmony.1500 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("FACS", "orig.ident", "sample"),
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
save(scobj.harmony.1500,
     scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
     file = "run2|nf1500_h_noreg scobj.harmony.Rdata")
  
# 2000 without regress cell cycle----
scobj.harmony.2000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("FACS", "orig.ident", "sample"),
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
save(scobj.harmony.2000,
     scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
     file = "run2|nf2000_h_noreg scobj.harmony.Rdata")

pdf("run2|nf1000_h_noreg clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
dev.off()

pdf("run2|nf1500_h_noreg clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
dev.off()

pdf("run2|nf2000_h_noreg clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
dev.off()
write.csv(data_frame(Finished = "TRUE"), "run2_noreg_fin.csv")
# regress cell cycle============================================================
# 1000 regccmt----
set.seed(1011)
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
library(harmony)
library(stringr)

load("~/escc_geo/2.downstream/run1|scobj4run2.Rdata")
scobj.harmony.1000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("FACS", "orig.ident", "sample"),
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
save(scobj.harmony.1000.20, file = "run2|nf1000_h_pc20_regccmt scobj.Rdata")
save(scobj.harmony.1000.25, file = "run2|nf1000_h_pc25_regccmt scobj.Rdata")
save(scobj.harmony.1000.30, file = "run2|nf1000_h_pc30_regccmt scobj.Rdata")
print("--->Finished!<---")

# 1500 regccmt----
set.seed(1011)
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
library(harmony)
library(stringr)

load("~/escc_geo/2.downstream/run1|scobj4run2.Rdata")
scobj.harmony.1500 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("FACS", "orig.ident", "sample"),
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
save(scobj.harmony.1500.20, file = "run2|nf1500_h_pc20_regccmt scobj.Rdata")
save(scobj.harmony.1500.25, file = "run2|nf1500_h_pc25_regccmt scobj.Rdata")
save(scobj.harmony.1500.30, file = "run2|nf1500_h_pc30_regccmt scobj.Rdata")
print("--->Finished!<---")

# 2000 regccmt----
set.seed(1011)
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
library(harmony)
library(stringr)

load("~/escc_geo/2.downstream/run1|scobj4run2.Rdata")
scobj.harmony.2000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("FACS", "orig.ident", "sample"),
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
save(scobj.harmony.2000.20, file = "run2|nf2000_h_pc20_regccmt scobj.Rdata")
save(scobj.harmony.2000.25, file = "run2|nf2000_h_pc25_regccmt scobj.Rdata")
save(scobj.harmony.2000.30, file = "run2|nf2000_h_pc30_regccmt scobj.Rdata")
print("--->Finished!<---")

# Run Umap ----
# noreg 1000 30 0.5
load("run2|nf1000_h_noreg scobj.harmony.Rdata")
scobj.harmony <- scobj.harmony.1000.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Visualization=================================================================
# scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.5")
pdf("run2|nf1000_h_pc30_noreg_res0.5 umap.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T) + NoLegend() 
DimPlot(scobj.harmony, group.by = "sample", reduction = "umap", label = T) + NoLegend() 
DimPlot(scobj.harmony, group.by = "FACS", reduction = "umap", label = T) + NoLegend() 
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.5", reduction = "umap", label = T)
dev.off()
# Annotation====================================================================
## Findmarkers------------------------------------------------------------------
Clist <- list()
n <- 1
for (i in 0:23) {
  Clist[[n]] <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.5", ident.1 = as.character(i))
  write.csv(Clist[[n]], paste0("C",n-1,".csv"))
  
  n <- n + 1
}

epithelial  <- c("EPCAM", "SFN", "KRT5", "KRT14") #, "SPRR3"
mesothelial <- c("PDPN", "MSLN", "UPK3B")
endothelial <- c("VWF", "PECAM1", "ENG", "CDH5") #, "CCL14"
fibroblast  <- c("FN1", "DCN", "COL1A1", "COL1A2") #, "COL3A1", "COL6A1"
SMC         <- c("TAGLN", "CNN1", "PRKG1", "FOXP2")
pericyte    <- c("RGS5", "ACTA2", "MYH11") # "MCAM", 
T_cell      <- c("CD2", "CD3D", "CD3E", "CD3G") 
B_cell      <- c("CD19", "CD79A", "MS4A1", "CD79B")
plasma      <- c("JCHAIN", "MZB1", "IGHG1", "SDC1") # , "CD79A"
myeloid     <- c("CD68", "CD163", "LYZ", "CD14","FCGR3A","C1QA","C1QB","S100A8","S100A9") # , "IL3RA","CST3", "LAMP3", "CLEC4C", "GCA"
Neutrophil  <- c("CSF3R","CXCL8","G0S2") # ,"IFITM2","FCGR3B","FPR1","BASP1","CXCR1","CXCR2","S100A11"
DC          <- c('CCR7',	'CLEC9A', 'CD1C',	'IRF7',	'LILRA4')
Mast        <- c('TPSAB1','CPA3','HPGDS','KIT') # 'VWA5A','SLC18A2','HDC','CAPG','RGS13','IL1RL1','FOSB','GATA2'
Cellcycle   <- c("PCNA","MKI67","DLGAP5","TOP2A")
schwann<- c("SOX10", "S100B")
neural <- c("TAC1", "CHAT",'PLP1', 'NRXN1', 'NRXN2', 'NRXN3')
arteries<- c("HEY1", "IGFBP3")
capillaries <- c("CD36", "CA4")
veins <- c("ACKR1")

p <- DotPlot(scobj.harmony, features = c(
  epithelial, mesothelial, endothelial, arteries, capillaries, veins, fibroblast, SMC, pericyte,
  T_cell, B_cell,plasma, myeloid, Neutrophil, Mast, DC, Cellcycle,
  schwann, neural
) %>% unique, 
group.by = "RNA_snn_res.0.5") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = "run2|marker.pdf", width=12, height=10, units="in")

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(0)] <- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(1)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(2)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(3)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(4)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(5)] <- "B cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(6)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(7)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(8)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(9)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(10)]<- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(11)]<- "B cell" # Plasma
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(12)]<- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(13)]<- "B cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(14)]<- "Mast"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(15)]<- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(16)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(17)]<- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(18)]<- "Myeloid" # DC
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(19)]<- "Myeloid" # DC
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(20)]<- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(21)]<- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(22)]<- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(23)]<- "Mesenchyme"

table(scobj.harmony$RNA_snn_res.0.5)
# 0     1     2     3     4     5     6     7     8     9    10    11    12 
# 24213 16087 14527 13493 12284 10983 10815  8809  8424  8393  4378  4251  2834 
# 13    14    15    16    17    18    19    20    21    22    23 
# 1611  1592  1314  1280   945   663   563   435   349   185   143 
table(scobj.harmony$cell_type)
# B cell Endothelial  Epithelial        Mast  Mesenchyme     Myeloid 
# 16845        9429       20519        1592       29798        9619 
# T cell 
# 60769

table(scobj.harmony$RNA_snn_res.0.5, scobj.harmony$FACS)
#     CD45_N CD45_P
# 0   24213      0
# 1       2  16085
# 2       0  14527
# 3       1  13492
# 4      21  12263
# 5       0  10983
# 6   10750     65
# 7    8809      0
# 8    4608   3816
# 9      10   8383
# 10      0   4378
# 11      3   4248
# 12   2834      0
# 13      1   1610
# 14      1   1591
# 15   1314      0
# 16   1264     16
# 17    940      5
# 18      0    663
# 19      1    562
# 20    435      0
# 21    349      0
# 22    185      0
# 23    143      0

meta      <- readxl::read_excel("../GSE160269_metadata.xlsx")
meta.data <- scobj.harmony@meta.data
meta.data <- meta.data %>% 
  left_join(meta, by = "sample")
rownames(meta.data) <- rownames(scobj.harmony@meta.data)
meta.data -> scobj.harmony@meta.data

# save(scobj.harmony, file = "run2|scobj.harmony(singlet_nf1000_h_noreg_pc30_res0.5)anno.Rdata")
load("run2|scobj.harmony(singlet_nf1000_h_noreg_pc30_res0.5)anno.Rdata")
