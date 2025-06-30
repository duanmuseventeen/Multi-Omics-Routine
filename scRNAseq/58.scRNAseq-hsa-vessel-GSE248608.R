rm(list = ls());gc()

# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
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
library(RColorBrewer)
library(slingshot) # Trajectory
# library(tradeSeq) # Trajectory
# library(monocle) # Trajectory
library(monocle3) # Trajectory
library("CellChat") # Communication
# library("iTALK") # Communication
# library(infercnv) # CNV
# library(copykat) # CNV
# library(ArchR)
# Load data---------------------------------------------------------------------
# options(future.globals.maxSize = 1000 * 1024^2) # 设置为1GB
set.seed(1011)

# https://github.com/friedpine/scRNASeq_GSE248608/blob/main/3c.R
mycol <-c("#DEEAB2","#64B473","#2D553C","#A1D8C9","#487C76","#7AAF93","#D0E4C6",
          "#F3746C","#BB4B95","#F8BFAF","#F7DDD3","#66CDF6","#598198","#D5E7F7",
          "#69B3CE","#D6D5EB","#7B8CBC","#7674AE", "#E3CEE4",
          "#AFB2B6","#C9BDB2","#F5C785","#ECAECF","#E9A943","#CAA57D","#A79388",
          "#EACF68","#F6F3A7","#C45337","#86382B","#EADCE4",
          "#EE5276","#9E6CAB","#74507B")

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

CA1 37487 4992 4840390
CA2 36601 4992 4936830
CA3 36601 4992 751821
NC1 37487 4992 5031030
NC2 36601 4992 2321360
NC3 36601 4992 1329964

# A single-cell suspension was prepared from the tibial anterior vessels isolated from 4 diabetic patients. 

hub <- c("CA1", "CA2", "CA3", "NC1", "NC2", "NC3")
names(hub) <- hub

setwd("GSE248608/")
CA1 <- "CA1" %>% Read10X %>% CreateSeuratObject(project = "CA1",assay = "RNA", min.cells = 0, min.features = 0)
CA2 <- "CA2" %>% Read10X %>% CreateSeuratObject(project = "CA2",assay = "RNA", min.cells = 0, min.features = 0)
CA3 <- "CA3" %>% Read10X %>% CreateSeuratObject(project = "CA3",assay = "RNA", min.cells = 0, min.features = 0)
NC1 <- "NC1" %>% Read10X %>% CreateSeuratObject(project = "NC1",assay = "RNA", min.cells = 0, min.features = 0)
NC2 <- "NC2" %>% Read10X %>% CreateSeuratObject(project = "NC2",assay = "RNA", min.cells = 0, min.features = 0)
NC3 <- "NC3" %>% Read10X %>% CreateSeuratObject(project = "NC3",assay = "RNA", min.cells = 0, min.features = 0)

all(intersect(rownames(CA1), rownames(CA2)) == intersect(rownames(NC1), rownames(NC2)))
# [1] TRUE

gene_symbol <- intersect(rownames(CA1), rownames(CA2))
length(gene_symbol)
# [1] 23705

sclist <- list(CA1, CA2, CA3, NC1, NC2, NC3)
GSE248608 <- merge(x = sclist[[1]], y = sclist[-1], add.cell.ids = hub)
# An object of class Seurat 
# 50383 features across 29952 samples within 1 assay 
# Active assay: RNA (50383 features, 0 variable features)
# 6 layers present: counts.CA1, counts.CA2, counts.CA3, counts.NC1, counts.NC2, counts.NC3

scobj <- GSE248608[gene_symbol,]
# An object of class Seurat 
# 23705 features across 29952 samples within 1 assay 
# Active assay: RNA (23705 features, 0 variable features)
# 6 layers present: counts.CA1, counts.CA2, counts.CA3, counts.NC1, counts.NC2, counts.NC3

scobj.rmENSG <- scobj[!str_detect(rownames(scobj), "^ENSG"),]
# An object of class Seurat 
# 23705 features across 29952 samples within 1 assay 
# Active assay: RNA (23705 features, 0 variable features)
# 6 layers present: counts.CA1, counts.CA2, counts.CA3, counts.NC1, counts.NC2, counts.NC3

dim(scobj.rmENSG@assays$RNA$counts)
# [1] 23705  4992

scobj <- scobj.rmENSG

table(scobj@meta.data$orig.ident)
# CA1  CA2  CA3  NC1  NC2  NC3 
# 4992 4992 4992 4992 4992 4992 
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(scobj, features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, features = c("percent.rp"), ncol = 1) + scale_y_continuous(breaks = c(10,20,25))
VlnPlot(scobj, features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nFeature = scobj@meta.data$nFeature_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nFeature, color = group)) +
  geom_density() +
  geom_vline(xintercept = c(200,300,400,500,1000,5000,7000,8000), color = "gray50", linetype = 2) +
  # scale_x_continuous(breaks = c(200,300,400,500,750,1000,3000,4000,5000)) +
  theme_classic()

data.frame(
  nCount_RNA = scobj@meta.data$nCount_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nCount_RNA, color = group)) +
  geom_density() +
  geom_vline(xintercept = c(0,2000,5000,10000,20000), color = "gray50", linetype = 2) +
  scale_x_continuous(limits = c(0,5000)) +
  theme_classic()

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 200 &
    # nFeature_RNA < 2500 & 
    # nCount_RNA < 2000 &
    percent.mt < 5
  # percent.rp < 20 &
  # percent.hb < 1
)

# An object of class Seurat 
# 23705 features across 21507 samples within 1 assay 
# Active assay: RNA (23705 features, 0 variable features)
# 6 layers present: counts.CA1, counts.CA2, counts.CA3, counts.NC1, counts.NC2, counts.NC3

table(scobj@meta.data$orig.ident)
# CA1  CA2  CA3  NC1  NC2  NC3 
# 3257 3202   13 3228 1894    7 

# Removing NC3 CA3 Load data----------------------------------------------------
CA1 <- "CA1" %>% Read10X %>% CreateSeuratObject(project = "CA1",assay = "RNA", min.cells = 0, min.features = 0)
CA2 <- "CA2" %>% Read10X %>% CreateSeuratObject(project = "CA2",assay = "RNA", min.cells = 0, min.features = 0)
NC1 <- "NC1" %>% Read10X %>% CreateSeuratObject(project = "NC1",assay = "RNA", min.cells = 0, min.features = 0)
NC2 <- "NC2" %>% Read10X %>% CreateSeuratObject(project = "NC2",assay = "RNA", min.cells = 0, min.features = 0)

all(intersect(rownames(CA1), rownames(CA2)) == intersect(rownames(NC1), rownames(NC2)))
# [1] TRUE

gene_symbol <- intersect(rownames(CA1), rownames(CA2))
length(gene_symbol)
# [1] 23705

sclist <- list(CA1, CA2, NC1, NC2)
GSE248608 <- merge(x = sclist[[1]], y = sclist[-1], add.cell.ids = hub)
# An object of class Seurat 
# 50383 features across 29952 samples within 1 assay 
# Active assay: RNA (50383 features, 0 variable features)
# 6 layers present: counts.CA1, counts.CA2, counts.CA3, counts.NC1, counts.NC2, counts.NC3

scobj <- GSE248608[gene_symbol,]
scobj.rmENSG <- scobj[!str_detect(rownames(scobj), "^ENSG"),]

scobj <- scobj.rmENSG

table(scobj@meta.data$orig.ident)
# CA1  CA2  NC1  NC2 
# 4992 4992 4992 4992
# Removing NC3 CA3 QC-----------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(scobj, features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, features = c("percent.rp"), ncol = 1) + scale_y_continuous(breaks = c(10,20,25))
VlnPlot(scobj, features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nFeature = scobj@meta.data$nFeature_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nFeature, color = group)) +
  geom_density() +
  geom_vline(xintercept = c(200,300,400,500,1000,5000,7000,8000), color = "gray50", linetype = 2) +
  # scale_x_continuous(breaks = c(200,300,400,500,750,1000,3000,4000,5000)) +
  theme_classic()

data.frame(
  nCount_RNA = scobj@meta.data$nCount_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nCount_RNA, color = group)) +
  geom_density() +
  geom_vline(xintercept = c(0,2000,5000,10000,20000), color = "gray50", linetype = 2) +
  scale_x_continuous(limits = c(0,10000)) +
  theme_classic()

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    nCount_RNA < 5000 &
    percent.mt < 10 &
    percent.rp < 25 &
    percent.hb < 1
)

# An object of class Seurat 
# 23705 features across 21507 samples within 1 assay 
# Active assay: RNA (23705 features, 0 variable features)
# 6 layers present: counts.CA1, counts.CA2, counts.CA3, counts.NC1, counts.NC2, counts.NC3

table(scobj@meta.data$orig.ident)
# CA1  CA2  NC1  NC2 
# 2351 3418 3336 3077
# Removing NC3 CA3 blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(GSE248608))
# [1] 246
scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj)
# [1] 23468 12182

any(rownames(scobj) %>% str_detect("^MT-")) # [1] FALSE
any(rownames(scobj) %>% str_detect("^RP[SL]")) # [1] FALSE

# save(scobj, file = "scobj.Rdata")
# Removing NC3 CA3 Nomalization & harmony & cluster----------------------------------------------
# Processes including cell cycle and tissue dissociation may influence the expression of not only the associated signature genes, but also the whole transcriptome (51).
# To minimize the impact of those processes, the S phase score and G2M phase score were calculated with function CellCycleScoring, and the DIG score was calculated with function AddModuleScore. 
# Then those scores were regressed out in the Seurat pipeline.

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

scobj.harmony.10 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.15 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.20 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.25 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.30 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

pdf("vst.2000.h clustree 10-30.pdf")
ElbowPlot(scobj.harmony, reduction = "pca", ndims = 50)
clustree(scobj.harmony.10, prefix = "RNA_snn_res.")
clustree(scobj.harmony.15, prefix = "RNA_snn_res.")
clustree(scobj.harmony.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.30, prefix = "RNA_snn_res.")
dev.off()

scobj.harmony <- scobj.harmony.15 %>% 
  RunUMAP(reduction = "harmony", dims = 1:15)
# Removing NC3 CA3 Visualization------------------------------------------------
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.2")
FeaturePlot(scobj.harmony, features = "contamination")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.2", reduction = "umap", label = T)
# Conclusion--------------------------------------------------------------------
# The sample NC3 and CA3 are low-qualitied. After removing them, the four samples
# cannot be clustered different clusters in UMAP.

# The Method described in the article is not reproducible.
# We filtered cells that have unique feature counts less than 200 or have > 5% 
# mitochondrial counts. Use scDblFinder to evaluate the doublet identification 
# cells [17], then perform normalizing and Scaling on the data. We chose 15
# Principal Component Analysis (PCA) dimensionality to find neighbors and clustered 
# data with 0.8 resolution.

