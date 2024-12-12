# Single-cell transcriptomics and surface epitope detection in human brain epileptic lesions identifies pro-inflammatory signaling
# GSE201048
# PMID: 35739273

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
scobj <- readRDS("C:/D/R project/Multi-Omics-Routine/scRNAseq/GSE201048_brain_imm_all_cells_25_cluster.rds")
scobj@assays$RNA$counts %>% dim
[1] 22950 85000
scobj@assays$ADT$counts %>% dim
[1]    16 85000
head(scobj@assays$ADT$counts[,1:6],16)

DimPlot(scobj, reduction = "tsne",
        group.by = "seurat_clusters",
        label = TRUE, pt.size = 0.5)
FeaturePlot(sce.all.annot, features = c("FAF2"))
VlnPlot(object = sce.all.annot, features = 'FAF2', slot = "counts", log = TRUE, group.by = "cell_type")

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
# QC---------------------------------------------------------
sce[["percent.mt"]] <- PercentageFeatureSet(sce, assay = "RNA", pattern = "^MT-")
sce[["percent.rp"]] <- PercentageFeatureSet(sce, assay = "RNA", pattern = "^RP[SL]")

pdf("P2-QC.pdf")
VlnPlot(sce, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(sce, features = c("nCount_RNA"), ncol = 1)
VlnPlot(sce, features = c("percent.mt"), ncol = 1)
VlnPlot(sce, features = c("percent.rp"), ncol = 1)
FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nFeature = sce@meta.data$nFeature_RNA
) %>% 
  ggplot(aes(x = nFeature)) +
  geom_density() + 
  scale_x_continuous(breaks = c(200,500,1000,3000,4000,5000))
dev.off()

# Cells with between 300 and 5,000 genes and mitochondrial percent reads less than 20 were kept for analysis.
sce <- subset(
  sce, 
  subset = nFeature_RNA > 300 & 
    nFeature_RNA < 5000 & 
    # nCount_RNA < 5000 &
    percent.mt < 20)  

pdf("P2-QC-2.pdf")
VlnPlot(sce, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(sce, features = c("nCount_RNA"), ncol = 1)
VlnPlot(sce, features = c("percent.mt"), ncol = 1)
VlnPlot(sce, features = c("percent.rp"), ncol = 1)
dev.off()

DefaultAssay(sce) <- 'RNA'
sce <- sce %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50)

# DefaultAssay(sce) <- 'ADT'
# sce <- sce %>% 
#   NormalizeData(normalization.method = 'CLR', margin = 2) %>% 
#   ScaleData() %>% 
#   RunPCA(features = rownames(sce), reduction.name = 'apca')

#Examine and visualize PCA results a few different ways
pdf("QC-PCA.pdf")
DimHeatmap(sce, dims = 1:6, cells = 500, balanced = TRUE)
ElbowPlot(sce, reduction = "pca", ndims = 50)
dev.off()
# Not Integrate Data-----------------------------------------
if(TRUE){
  # Principal-component (PC) analysis was performed   on the 2,000 most variable genes, and the first 20 PCs were used for t-SNE and UMAP for data embedding into two dimensions.
  
  # Cluster---------------------------------------------------------------------
  pdf("clustree without integrating data.pdf")
  sce <- FindNeighbors(sce, reduction = "pca", dims = 1:20)
  sce <- FindClusters(sce, resolution = seq(0.1, 1, 0.1))
  clustree(sce, prefix = "RNA_snn_res.")
  dev.off()
  
  # Visualization---------------------------------------------------------------
  #UMAP
  sce <- RunUMAP(sce, reduction = "pca", min_dist = 0.3, dims = 1:20)
  #T-SNE
  sce <- RunTSNE(sce, reduction = "pca", dims = 1:20)
  
  pdf("umap and tsne without integrating data.pdf")
  sce <- RegroupIdents(sce, metadata = "RNA_snn_res.0.9")
  DimPlot(sce, group.by = "orig.ident", reduction = "umap", label = T)
  DimPlot(sce, group.by = "RNA_snn_res.0.9", reduction = "umap", label = T)
  DimPlot(sce, group.by = "orig.ident", reduction = "tsne", label = T)
  DimPlot(sce, group.by = "RNA_snn_res.0.9", reduction = "tsne", label = T)
  dev.off()
}
# Annotation--------------------------------------------------------------------
# wsnn_res.0.9
sce <- RegroupIdents(sce, metadata = "wsnn_res.0.6")
sce.copy <- sce

# Cell clustering into subsets using graph-based Louvain clustering
# algorithms resulted in 26 clusters. Based on CD45 expression
# levels, 13 clusters (0–7, 9–12 and 14) were identified as microglia
# (CD45lo) and six clusters (8, 15–17, 19 and 21) were identified as
# infiltrating immune cells (CD45hi). All other clusters were CD45−
# and were identified using marker gene expression levels. Clusters
# 13, 20 and 22 expressed genes (CLDN5, MYH11, ABCC9, VWF and
# ACTA2) specific to cells of the neurovascular unit (NVU)15, while
# cluster 18 expressed genes (MAG and MOG) specific to oligodendrocytes
# (Supplementary Fig. 1). Among the immune cell subsets,
# cluster 19 had B cell marker proteins (CD19 and CD20) and cluster
# 17 had macrophage marker proteins (CD45hi, CD14 and CD11b).
# Clusters 8, 15 and 16 were identified as T cell clusters (CD45hiCD3+)
# (Fig. 1a). Infiltrating immune cells were observed in all 11 tissues
# from the olfactory, frontal or temporal lobes, irrespective of
# their location in the brain (Fig. 1c).

Lymphoid <- c("CD45")
microglial <- c("CX3CR1", "P2RY12", "TREM2")
neurovascular <- c("CLDN5", "MYH11", "ABCC9", "VWF", "ACTA2")
oligodendrocytes <- c("MAG", "MOG")
B <- c("CD19", "CD20")
M <- ("CD45hi", "CD14", "CD11b")
Tcell <- c("CD45", "CD3")

if(TRUE){
  # pdf("Vln.pdf")
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = ACM, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = VCM, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = EC, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = FB, layer = "counts", log = TRUE)
  # dev.off()
}

# sce@meta.data$cell_type <- "Unknown" #
# sce@meta.data$cell_type[sce@meta.data$wsnn_res.0.6 == 0] <- ""
# 
# p1 <- DimPlot(sce, reduction = "umap",
#               group.by = "wsnn_res.0.6",
#               label = TRUE, pt.size = 0.5) 
# p2 <- DimPlot(sce, reduction = "umap",
#               group.by = "cell_type",
#               label = TRUE, pt.size = 0.5)
# p3 <- DimPlot(sce, reduction = "umap",
#               group.by = "orig.ident",
#               label = TRUE, pt.size = 0.5) 
# p1 + p2 + p3
# # ggsave(..., width = 18, height = 6)
# 
# save(sce, file = "C:/D/.../sce.annot.Rdata")


