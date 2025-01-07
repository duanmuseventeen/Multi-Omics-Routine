# Single-cell transcriptomics and surface epitope detection in human brain epileptic lesions identifies pro-inflammatory signaling
# GSE201048
# PMID: 35739273

# https://satijalab.org/seurat/articles/seurat5_integration.html#:~:text=Seurat%20v5%20assays%20store%20data%20in%20layers.%20These,%28layer%3D%27counts%27%29%2C%20normalized%20data%20%28layer%3D%27data%27%29%2C%20or%20z-scored%2Fvariance-stabilized%20data%20%28layer%3D%27scale.data%27%29.
# Note that since the data is split into layers, normalization and variable feature identification is performed for each batch independently (a consensus set of variable features is automatically identified).
# Once integrative analysis is complete, you can rejoin the layers - which collapses the individual datasets together and recreates the original counts and data layers. You will need to do this before performing any differential expression analysis. However, you can always resplit the layers in case you would like to reperform integrative analysis.

# https://epicimmuneatlas.org/NatNeu2022/
# https://zenodo.org/records/6477100

# Reference---------------------------------------------------------------------
# https://zhuanlan.zhihu.com/p/567253121
# https://www.jianshu.com/p/442d890d12ea # about cite-seq
# https://cite-seq.com/ # about cite-seq
# https://www.bilibili.com/video/BV1xA411v7MR/?spm_id_from=333.999.0.0&vd_source=7687d926bc3779b048e880123e66bca4
# https://www.jianshu.com/p/0a28b293d6dd
# An entropy-based metricfor assessing the purity of single cell populations
# https://github.com/PaulingLiu/ROGUE
# https://zhuanlan.zhihu.com/p/667102765
# https://satijalab.org/seurat/articles/integration_introduction
# https://www.jianshu.com/p/7c43dc99c4b1
# https://blog.csdn.net/qq_52813185/article/details/134738735 # graphics
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
# https://cloud.tencent.com/developer/article/1865543
# https://www.jianshu.com/p/d1394ee923a2
# https://www.jianshu.com/p/c438d545f696

rm(list = ls()); gc()
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree) # use for determine the optimal resolution
library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
library(clustree)
library(ROGUE)
library(scDblFinder)
library(CellMixS)
# Observe-------------------------------------------------------------------------
scobj <- readRDS("GSE201048_brain_imm_all_cells_25_cluster.rds")

scobj@assays$RNA$counts %>% dim
# [1] 22950 85000
scobj@assays$RNA$data %>% dim
# [1] 22950 85000
scobj@assays$RNA$scale.data %>% dim
# [1]    50 85000

scobj@assays$ADT$counts %>% dim
# [1]    16 85000
scobj@assays$ADT$data %>% dim
# [1]    16 85000
scobj@assays$ADT$scale.data %>% dim
# [1]    16 85000

head(scobj@assays$ADT$counts[,1:6],16)

table(scobj@meta.data$sampleName)
# P5.B  P5.A  P3.A    P2  P1.B  P1.A    P4  P3.B  P6.A  P6.B  P6.C 
# 7074 10281  7783  6842  4464  6429  4596 11153  9099  8496  8783

DimPlot(scobj, reduction = "tsne",
        group.by = "seurat_clusters",
        label = TRUE, pt.size = 0.5)
FeaturePlot(sce.all.annot, features = c("FAF2"))
VlnPlot(object = sce.all.annot, features = 'FAF2', slot = "counts", log = TRUE, group.by = "cell_type")

if(TRUE){
  scobj@commands
  
  # NormalizeData(pbmc)
  # FindVariableFeatures(pbmc, nfeatures = 2000)
  # ScaleData(pbmc)
  # RunPCA(pbmc, verbose = FALSE)
  # FindNeighbors(pbmc, dims = 1:20)
  # FindClusters(pbmc, resolution = 0.8)
  # RunUMAP(pbmc, dims = 1:20)
  # RunTSNE(pbmc, dims = 1:20, method = "FIt-SNE")
  # 
  # NormalizeData(pbmc_protein, assay = "ADT", normalization.method = "CLR")
  # ScaleData(pbmc_protein, assay = "ADT")
}

# Linux-------------------------------------------------------------------------
for i in $(seq 1 11);
do
mkdir Sample${i}
mv GSM*_Sample${i}_* Sample${i}
done

for i in $(seq 1 11);
do
mv Sample${i}/*features.tsv.gz Sample${i}/features.tsv.gz
mv Sample${i}/*barcodes.tsv.gz Sample${i}/barcodes.tsv.gz
mv Sample${i}/*matrix.mtx.gz Sample${i}/matrix.mtx.gz
done

ls -lh Sample*/

# To match rows between samples
# cellranger输出的matrix.mtx除了%以外的第一行的三个数分别代表：基因、barcode、matrix.mtx矩阵的行数
for i in $(seq 1 11);
  do
  gunzip Sample${i}/features.tsv.gz 
  sed -i -e "/P298/d" -i -e "/N212/d" Sample${i}/features.tsv
  gzip Sample${i}/features.tsv
  
  gunzip Sample${i}/barcodes.tsv.gz
  sed -i -e "/32755/d" -i -e "/32756/d" Sample${i}/barcodes.tsv
  gzip Sample${i}/barcodes.tsv
  done
  
for i in $(seq 1 11);
  do
  gunzip Sample${i}/matrix.mtx.gz 
  sed -i -e "3s/32756/32754/" Sample${i}/matrix.mtx
  sed -i -e "/32755/d" -i -e "/32756/d" Sample${i}/matrix.mtx
  gzip Sample${i}/matrix.mtx
  done
# Load data----------------------------------------------------------------------
# # sample5 sample6 have 16 ADTs, whereas others have 18 ADTs
# rownames(sample1[[2]])
# [1] "CD3_TotalSeqB"    "CD4_TotalSeqB"    "CD8a_TotalSeqB"   "CD14_TotalSeqB"
# [5] "CD16_TotalSeqB"   "CD19_TotalSeqB"   "CD25_TotalSeqB"   "CD56_TotalSeqB"
# [9] "CD20_TotalSeqB"   "CD335_TotalSeqB"  "CD69_TotalSeqB"   "CD197_TotalSeqB"
# [13] "CD27_TotalSeqB"   "HLA-DR_TotalSeqB" "CD11b_TotalSeqB"  "CD45_TotalSeqB"
# [17] "prip_TotalSeqB"   "n212_TotalSeqB"
# 
# rownames(sample5[[2]])
# [1] "CD3_TotalSeqB"    "CD4_TotalSeqB"    "CD8a_TotalSeqB"   "CD14_TotalSeqB"
# [5] "CD16_TotalSeqB"   "CD19_TotalSeqB"   "CD25_TotalSeqB"   "CD56_TotalSeqB"
# [9] "CD20_TotalSeqB"   "CD335_TotalSeqB"  "CD69_TotalSeqB"   "CD197_TotalSeqB"
# [13] "CD27_TotalSeqB"   "HLA-DR_TotalSeqB" "CD11b_TotalSeqB"  "CD45_TotalSeqB"

hub <- paste0("Sample", c(1:11))
names(hub) <- hub

if(TRUE){
  # match rows in R
  sceList.GEX <- lapply(hub, function(x){
    sc10X <- Read10X(x)
    sce <- CreateSeuratObject(
      counts = sc10X$`Gene Expression`,
      project = x,
      min.cells = 0,
      min.features = 0,
      assay = "RNA")
    return(sce)
  })
  
  sceList.ADT <- lapply(hub, function(x){
    sc10X <- Read10X(x)
    sce <- CreateSeuratObject(
      counts = sc10X$`Antibody Capture`[1:16,],
      project = x,
      min.cells = 0,
      min.features = 0,
      assay = "ADT")
    return(sce)
  })
  
  scobj.GEX <- merge(x = sceList.GEX[[1]], y = sceList.GEX[-1], add.cell.ids = hub)
  scobj.ADT <- merge(x = sceList.ADT[[1]], y = sceList.ADT[-1], add.cell.ids = hub)
  
  table(scobj.GEX@meta.data$orig.ident)
  table(scobj.ADT@meta.data$orig.ident)
  Sample1 Sample10 Sample11  Sample2  Sample3  Sample4  Sample5  Sample6
  6777     8986     9164     5148     7346     8360    11730     5197
  Sample7  Sample8  Sample9
  10601     7541     9491
  
  Sample1 Sample10 Sample11  Sample2  Sample3  Sample4  Sample5  Sample6
  6777     8986     9164     5148     7346     8360    11730     5197
  Sample7  Sample8  Sample9
  10601     7541     9491
  
  scobj.GEX[["ADT"]] <- scobj.ADT[["ADT"]]
  
  sce <- scobj.GEX
}

if(FALSE){
  # match rows in shell
  sc10X <- Read10X(hub)
  
  sce <- CreateSeuratObject(
    counts = sc10X[["Gene Expression"]],
    min.cells = 0,
    min.features = 0,
    project = "GEX",
    assay = "RNA")
  
  adt <- CreateAssayObject(sc10X[["Antibody Capture"]])
  all.equal(colnames(sc10X[["Gene Expression"]]), colnames(sc10X[["Antibody Capture"]]))
  # [1] TRUE
  sce[["ADT"]] <- adt
}

table(sce@meta.data$orig.ident)
# Add meta.data-----------------------------------------------------------------
# QC----------------------------------------------------------------------------
sce[["percent.mt"]] <- PercentageFeatureSet(sce, assay = "RNA", pattern = "^MT-")
sce[["percent.rp"]] <- PercentageFeatureSet(sce, assay = "RNA", pattern = "^RP[SL]")
sce[["percent.hb"]] <- PercentageFeatureSet(sce, pattern = "^HB[^(P)]")

# pdf("P2-QC.pdf")
# VlnPlot(sce, features = c("nFeature_RNA"), ncol = 1)
# VlnPlot(sce, features = c("nCount_RNA"), ncol = 1)
# VlnPlot(sce, features = c("percent.mt"), ncol = 1)
# VlnPlot(sce, features = c("percent.rp"), ncol = 1)
# VlnPlot(sce, features = c("percent.hb"), ncol = 1)
# FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# data.frame(
#   nFeature = sce@meta.data$nFeature_RNA
# ) %>% 
#   ggplot(aes(x = nFeature)) +
#   geom_density() + 
#   scale_x_continuous(breaks = c(200,500,1000,3000,4000,5000))
# dev.off()

# Cells with between 300 and 5,000 genes and mitochondrial percent reads less than 20 were kept for analysis.
sce <- subset(
  sce, 
  subset = nFeature_RNA > 500 & 
    nFeature_RNA < 4500 & 
    nCount_RNA < 15000 &
    percent.mt < 10 &
    percent.hb < 1)  

# pdf("P2-QC-2.pdf")
# VlnPlot(sce, features = c("nFeature_RNA"), ncol = 1)
# VlnPlot(sce, features = c("nCount_RNA"), ncol = 1)
# VlnPlot(sce, features = c("percent.mt"), ncol = 1)
# VlnPlot(sce, features = c("percent.rp"), ncol = 1)
# VlnPlot(sce, features = c("percent.hb"), ncol = 1)
# FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# dev.off()

table(sce@meta.data$orig.ident)
# Sample1 Sample10 Sample11  Sample2  Sample3  Sample4  Sample5  Sample6
# 5693     7229     7813     3616     5856     6292     9509     3056
# Sample7  Sample8  Sample9
# 9117     6008     8172

# We first perform pre-processing and dimensional reduction on both assays independently. We use standard normalization, but you can also use SCTransform or any alternative method.
# https://satijalab.org/seurat/articles/seurat5_integration.html#:~:text=Seurat%20v5%20assays%20store%20data%20in%20layers.%20These,%28layer%3D%27counts%27%29%2C%20normalized%20data%20%28layer%3D%27data%27%29%2C%20or%20z-scored%2Fvariance-stabilized%20data%20%28layer%3D%27scale.data%27%29.
# Note that since the data is split into layers, normalization and variable feature identification is performed for each batch independently (a consensus set of variable features is automatically identified).
# Once integrative analysis is complete, you can rejoin the layers - which collapses the individual datasets together and recreates the original counts and data layers. You will need to do this before performing any differential expression analysis. However, you can always resplit the layers in case you would like to reperform integrative analysis.

DefaultAssay(sce) <- 'RNA'
sce <- sce %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50)

DefaultAssay(sce) <- 'ADT'
sce <- sce %>%
  NormalizeData(normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>%
  RunPCA(features = rownames(sce), npcs = 50, reduction.name = 'apca')
# 9 apca

pdf("PCA-samples.pdf")
DimPlot(sce.harmony, group.by = "orig.ident", reduction = "pca", label = T)
dev.off()
# Integration-------------------------------------------------------------------
DefaultAssay(sce) <- 'RNA'

sce.harmony <- RunHarmony(
  sce,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  reduction.save = "harmony.pca") %>% 
  RunHarmony(
  group.by.vars = "orig.ident",
  reduction.use = "apca",
  reduction.save = "harmony.apca")

sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
sce[["ADT"]] <- JoinLayers(sce[["ADT"]])
# Identify multimodal neighbors.------------------------------------------------
# These will be stored in the neighbors slot, 
# and can be accessed using sce[['weighted.nn']]
# The WNN graph can be accessed at sce[["wknn"]], 
# and the SNN graph used for clustering at sce[["wsnn"]]
# Cell-specific modality weights can be accessed at sce$RNA.weight
sce.harmony <- FindMultiModalNeighbors(
  sce.harmony, 
  reduction.list = list("harmony.pca", "harmony.apca"), 
  dims.list = list(1:30, 1:9), 
  modality.weight.name = c("RNA.weight", "ADT.weight")) %>% 
# Cluster & Visualization-------------------------------------------------------
  FindClusters(
    graph.name = "wsnn", 
    algorithm = 3,
    resolution = seq(0.01,0.1,0.01), 
    verbose = FALSE)

pdf("clustree.pdf", width = 10, height = 12)
clustree(sce.harmony, prefix = "wsnn_res.")
dev.off()

# RunUMAP(
#   sce.harmony, 
#   dims = 1:30,
#   nn.name = "weighted.nn",
#   reduction = "harmony.pca",
#   reduction.name = "wnn.umap", 
#   reduction.key = "wnnUMAP_")
# Error in RunUMAP.Seurat(sce.harmony, dims = 1:30, nn.name = "weighted.nn",  :
#                           Only one parameter among 'dims', 'nn.name', 'graph', or 'features' should be used at a time to run UMAP
#                         
sce.harmony <- RunUMAP(
  sce.harmony, 
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap", 
  reduction.key = "wnnUMAP_")
sce.harmony <- RunTSNE(
  sce.harmony, 
  nn.name = "weighted.nn", 
  reduction.name = "wnn.tsne", 
  reduction.key = "wnnTSNE_")

pdf("clusters.pdf")
DimPlot(sce.harmony, group.by = "orig.ident", reduction = "wnn.umap", label = T)
DimPlot(sce.harmony, group.by = "wsnn_res.0.02", reduction = "wnn.umap", label = T)
# DimPlot(sce.harmony, group.by = "wsnn_res.0.02", reduction = "wnn.tsne", label = T)
dev.off()
# Cell cycle--------------------------------------------------------------------
exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sce.harmony <- CellCycleScoring(
  sce.harmony,
  s.features = s.genes, 
  g2m.features = g2m.genes, 
  set.ident = TRUE)

pdf("Cell cycle.pdf")
sce.harmony <- RunPCA(sce.harmony, features = c(s.genes, g2m.genes), reduction.name = "cell_cycle")
DimPlot(sce.harmony, reduction = "cell_cycle")
DimPlot(sce.harmony, reduction = "wnn.umap")
dev.off()
# Filter doublets---------------------------------------------------------------
if(TRUE){
  # scDblFinder
  sce.harmony <- as.SingleCellExperiment(sce.harmony, assay = "RNA")
  sce.harmony.scdbl <- scDblFinder(
    sce.harmony, 
    clusters = sce.harmony$wsnn_res.0.02,
    samples = sce.harmony$orig.ident,
    dbr.sd = 1) %>% 
    as.Seurat %>%
    subset(scDblFinder.class == "singlet")
  
  dim(sce.harmony)
  # [1] 32738 72361
  dim(sce.harmony.scdbl)
  # [1] 32738 66751
}
# a serious wf requires removing the doublets and re-log-normalize, integrate,----
# cluster, visualize-------------------------------------------------------------
sce <- scobj.GEX
# Add meta.data-----------------------------------------------------------------
# QC----------------------------------------------------------------------------
sce[["percent.mt"]] <- PercentageFeatureSet(sce, assay = "RNA", pattern = "^MT-")
sce[["percent.rp"]] <- PercentageFeatureSet(sce, assay = "RNA", pattern = "^RP[SL]")
sce[["percent.hb"]] <- PercentageFeatureSet(sce, pattern = "^HB[^(P)]")

sce <- subset(
  sce, 
  subset = nFeature_RNA > 500 & 
    nFeature_RNA < 4500 & 
    nCount_RNA < 15000 &
    percent.mt < 10 &
    percent.hb < 1)  

sce@meta.data$id <- rownames(sce@meta.data) 
sce@meta.data <- sce@meta.data %>% 
  mutate(singlet = case_when(id %in% rownames(sce.harmony.scdbl@meta.data) ~ 1,
                             TRUE ~ 0))

sce <- subset(sce, subset = singlet == 1)

DefaultAssay(sce) <- 'RNA'
sce <- sce %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50)

DefaultAssay(sce) <- 'ADT'
sce <- sce %>%
  NormalizeData(normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>%
  RunPCA(features = rownames(sce), npcs = 50, reduction.name = 'apca')
# 9 apca
# Integration-------------------------------------------------------------------
DefaultAssay(sce) <- 'RNA'

sce.harmony <- RunHarmony(
    sce,
    group.by.vars = "orig.ident",
    reduction.use = "pca",
    reduction.save = "harmony.pca") %>% 
  RunHarmony(
    group.by.vars = "orig.ident",
    reduction.use = "apca",
    reduction.save = "harmony.apca") 
sce.harmony[["RNA"]] <- JoinLayers(sce.harmony[["RNA"]])
sce.harmony[["ADT"]] <- JoinLayers(sce.harmony[["ADT"]])
# Cluster & Visualization-------------------------------------------------------
sce.harmony <- sce.harmony %>% 
  FindMultiModalNeighbors(
    reduction.list = list("harmony.pca", "harmony.apca"), 
    dims.list = list(1:30, 1:9), 
  modality.weight.name = c("RNA.weight", "ADT.weight")) %>%
    FindClusters(
    graph.name = "wsnn", 
    algorithm = 3,
    resolution = seq(0.01,0.1,0.01), 
    verbose = FALSE)

pdf("clustree-R2.pdf", width = 10, height = 12)
clustree(sce.harmony, prefix = "wsnn_res.")
dev.off()

sce.harmony <- RunUMAP(
  sce.harmony, 
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap", 
  reduction.key = "wnnUMAP_")
sce.harmony <- RunTSNE(
  sce.harmony, 
  nn.name = "weighted.nn", 
  reduction.name = "wnn.tsne", 
  reduction.key = "wnnTSNE_")

pdf("clusters.pdf")
DimPlot(sce.harmony, group.by = "orig.ident", reduction = "wnn.umap", label = T)
DimPlot(sce.harmony, group.by = "wsnn_res.0.05", reduction = "wnn.umap", label = T)
# DimPlot(sce.harmony, group.by = "wsnn_res.0.05", reduction = "wnn.tsne", label = T)
dev.off()

table(sce.harmony@meta.data$wsnn_res.0.05)
# 0     1    10    11    12    13    14    15    16    17    18    19     2    20 
# 21201 13910   714   442     4     3     3     2     2     2     2     2  7111     2 
# 21    22    23    24    25    26    27    28    29     3     4     5     6     7 
# 2     2     2     2     2   205     2     2     2  6796  3948  3340  3254  2897 
# 8     9 
# 1487  1408 

save(sce.harmony, file = "GSE201048_without_cell_type.Rdata")
