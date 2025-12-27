# Single-nucleus RNA sequencing in ischemic cardiomyopathy reveals common transcriptional profile underlying end-stage heart failure
# 36790929

# Reference---------------------------------------------------------------------
# https://www.jianshu.com/p/0a28b293d6dd
# An entropy-based metricfor assessing the purity of single cell populations
# https://github.com/PaulingLiu/ROGUE
# https://zhuanlan.zhihu.com/p/667102765
# https://satijalab.org/seurat/articles/integration_introduction
# https://www.jianshu.com/p/7c43dc99c4b1
# https://blog.csdn.net/qq_52813185/article/details/134738735 # graphics

# Load pkgs---------------------------------------------------------------------
rm(list = ls())

library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(clustree) # use for determine the optimal resolution
# library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)

# Load data---------------------------------------------------------------------
# Download from # https://singlecell.broadinstitute.org/single_cell/study/SCP1849/single-nucleus-rna-sequencing-in-ischemic-cardiomyopathy-reveals-common-transcriptional-profile-underlying-end-stage-heart-failure?genes=FSHR&cluster=ICM_UMAP_V1.txt&spatialGroups=--&annotation=biosample_id--group--study&subsample=all&tab=annotatedScatter#study-visualize
h.ICM.l.v <- Read10X(data.dir = "ICM/")
h.ICM.l.v <- CreateSeuratObject(counts = h.ICM.l.v, 
                                project = "Human_ICM_left_ventricle", 
                                min.cells = 3,
                                min.features = 200)
meta <- readxl::read_excel("ICM_MetaData_V1.xlsx")

meta <- data.frame(
  NAME = colnames(h.ICM.l.v),
  stringsAsFactors = FALSE) %>% 
  left_join(meta, by = "NAME")

all(meta$NAME == colnames(h.ICM.l.v))
# [1] TRUE

h.ICM.l.v@meta.data$orig.ident <- meta$donor_id
h.ICM.l.v@meta.data$barcodes <- meta$NAME
h.ICM.l.v@meta.data$sample <- meta$donor_id
h.ICM.l.v@meta.data$cell_type_pub <- meta$cell_type_leiden05
h.ICM.l.v@meta.data$sex <- meta$sex
h.ICM.l.v@meta.data$disease <- meta$disease__ontology_label

# QC----------------------------------------------------------------------------
h.ICM.l.v[["percent.mt"]] <- PercentageFeatureSet(h.ICM.l.v, pattern = "^MT-")
h.ICM.l.v[["percent.rp"]] <- PercentageFeatureSet(h.ICM.l.v, pattern = "^RP[SL]")
VlnPlot(h.ICM.l.v, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(h.ICM.l.v, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

h.ICM.l.v <- subset(h.ICM.l.v, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 5 & percent.rp < 1)  
VlnPlot(h.ICM.l.v, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)

h.ICM.l.v <- h.ICM.l.v %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50)

#Examine and visualize PCA results a few different ways
print(h.ICM.l.v[["pca"]], dims = 1:10, nfeatures = 5)
VizDimLoadings(h.ICM.l.v, dims = 1:2, reduction = "pca")
DimPlot(h.ICM.l.v, reduction = "pca")
DimHeatmap(h.ICM.l.v, dims = 1:4, cells = 500, balanced = TRUE)
# Perform analysis without integration (just try)------------------------------------------
if(FALSE){
  h.ICM.l.v <- FindNeighbors(h.ICM.l.v, dims = 1:50)
  h.ICM.l.v <- FindClusters(h.ICM.l.v, resolution = seq(0.1,0.8,by = 0.1)) 
  
  #UMAP
  h.ICM.l.v <- RunUMAP(h.ICM.l.v, dims = 1:50)
  # head(h.ICM.l.v@reductions$umap@cell.embeddings)
  #note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
  p1 <- DimPlot(h.ICM.l.v, group.by = "RNA_snn_res.0.1", reduction = "umap", label = T)
  #T-SNE
  h.ICM.l.v <- RunTSNE(h.ICM.l.v, dims = 1:50, check_duplicates = FALSE)
  # head(h.ICM.l.v@reductions$tsne@cell.embeddings)
  p2 <- DimPlot(h.ICM.l.v, reduction = "tsne", label = T)
  # p1 + p2
  
  # select optimal res
  # 1
  p1 <- clustree(h.ICM.l.v, prefix = "RNA_snn_res.")
  p2 <- DimPlot(h.ICM.l.v, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T)
  p1 + p2
  # 2
  expr <- h.ICM.l.v@assays$RNA$counts %>% matrix
  p1 <- rogue(expr, labels = h.ICM.l.v@meta.data$RNA_snn_res.0.1, 
              samples = meta$donor_id, platform='UMI', span =0.6) %>%
    rogue.boxplot()
  p2 <- rogue(expr, labels = h.ICM.l.v@meta.data$RNA_snn_res.0.3, 
              samples = meta$donor_id, platform='UMI', span=0.6) %>%
    rogue.boxplot()
  p1 + p2
  
  # save(h.ICM.l.v, file = "h.ICM.l.v.Rdata") 
  # 
  # load("h.ICM.l.v.Rdata")
  
  # RNA_snn_res.0.3
  ACM <- c('TNNT2', 'TTN', 'MYH6', 'MYH7') # 5 6 7 8 9 10 11 13 14 18
  VCM <- c('MYH7', 'MYL2', 'FHL2') # 5 6 7 8 9 10 11 13 14 18
  EC <- c('VWF', 'CDH5', 'THBD', 'TEK', 'ENG', 'PECAM1') # 3 12 16 ?15
  FB <- c('FN1', 'VIM', 'COL1A1', 'COL1A2') # 0 1 22 23
  M <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74') # 2 21
  Myeloid	<- c('CD163', 'CD14', 'C1QA') # 2 21
  P <- c('HCN1', 'HCN4', 'CACNA1D') # pacemaker cells (P cells)
  Peri <- c('RGS5', 'ABCC9', 'KCNJ8') # 4 
  SMC <- c('MYH11', 'ACTA2','CNN1', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') 
  Meso <- c('MSLN', 'WT1', 'BNC1') 
  Neural <- c('PLP1', 'NRXN1', 'NRXN3') # 20
  Adi <- c('GPAM', 'FASN', 'LEP') # 17
  Lymph <- c('CD3E', 'IL7R', 'CD40LG') # 19
  Mast <- c('KIT', 'CPA3')
  RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ')
  
  VlnPlot(h.ICM.l.v, group.by = "RNA_snn_res.0.3", features = EC, slot = "counts", log = TRUE)
  unknown.markers %>% arrange(desc(avg_log2FC)) %>% head(30)
  
  new.cluster.ids <- c("FB", "FB", "Myeloid", "EC", "Peri", "CM", 
                       "CM", "CM", "CM", "CM", "CM", 
                       "CM", "EC", "CM", "CM", "Unknown",
                       "Myeloid", "Adi", "CM", "Lymphoid", "Neural",
                       "Myeloid", "FB", "Meso")
  names(new.cluster.ids) <- levels(h.ICM.l.v)
  h.ICM.l.v <- RenameIdents(h.ICM.l.v, new.cluster.ids)
  DimPlot(h.ICM.l.v, reduction = "umap", 
          # group.by = "orig.ident", 
          label = TRUE, pt.size = 0.5) + NoLegend()
}
# Integration bt harmony--------------------------------------------------------
# Integration first, then FindNeighbors and FindClusters
# UMAP is used more commonly
sce.all.harmony <- RunHarmony(h.ICM.l.v, 
                              group.by.vars = "orig.ident", 
                              reduction.use = "pca",
                              reduction.save = "harmony")

sce.all.harmony <- FindNeighbors(sce.all.harmony, reduction = "harmony", dims = 1:50)
sce.all.harmony <- FindClusters(sce.all.harmony, resolution = seq(0.1, 0.8, 0.1))

#UMAP
sce.all.harmony <- RunUMAP(sce.all.harmony, reduction = "harmony", dims = 1:50)
p1 <- DimPlot(sce.all.harmony, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T)
#T-SNE
sce.all.harmony <- RunTSNE(sce.all.harmony, reduction = "harmony", dims = 1:50, check_duplicates = FALSE)
p2 <- DimPlot(sce.all.harmony, reduction = "tsne", label = T)
p1 + p2

DimPlot(sce.all.harmony, group.by = "sample", reduction = "umap", label = T)
DimPlot(sce.all.harmony, group.by = "sample", reduction = "tsne", label = T)

# Select optimal res------------------------------------------------------------
# 1
clustree(sce.all.harmony, prefix = "RNA_snn_res.")
# # 2
# expr <- sce.all.harmony@assays$RNA$counts %>% matrix
# p1 <- rogue(expr, labels = sce.all.harmony@meta.data$RNA_snn_res.0.1, 
#             samples = meta$donor_id, platform='UMI', span =0.6) %>%
#   rogue.boxplot()
# p2 <- rogue(expr, labels = sce.all.harmony@meta.data$RNA_snn_res.0.3, 
#             samples = meta$donor_id, platform='UMI', span=0.6) %>%
#   rogue.boxplot()
# p1 + p2

save(sce.all.harmony, file = "sce.all.harmony.Rdata")
# Annotate Cell Type------------------------------------------------------------
# RNA_snn_res.0.3
sce.all.harmony <- RegroupIdents(sce.all.harmony, metadata = "RNA_snn_res.0.3")
sce.all.harmony.copy <- sce.all.harmony

ACM <- c('TNNT2', 'TTN', 'MYH6', 'MYH7') # 1 2 ?11
VCM <- c('MYH7', 'MYL2', 'FHL2') # 1 2 ?11
EC <- c('VWF', 'CDH5', 'THBD', 'TEK', 'ENG', 'PECAM1') # 3 14
FB <- c('FN1', 'VIM', 'COL1A1', 'COL1A2') # 0 12
M <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74') # 4 13
Myeloid	<- c('CD163', 'CD14', 'C1QA') # 4 13
P <- c('HCN1', 'HCN4', 'CACNA1D') # pacemaker cells (P cells) 
Peri <- c('RGS5', 'ABCC9', 'KCNJ8') # ?1 ?2 5
SMC <- c('MYH11', 'ACTA2','CNN1', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') # 6
Meso <- c('MSLN', 'WT1', 'BNC1') # ?11
Neural <- c('PLP1', 'NRXN1', 'NRXN3') # 10
Adi <- c('GPAM', 'FASN', 'LEP') # 8
Lymph <- c('CD3E', 'IL7R', 'CD40LG') # 9
Mast <- c('KIT', 'CPA3')
NKT <- c('SLAMF1', 'SLAMF6', 'TGFBR', 'VA24', 'JA18') # 
RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ')

VlnPlot(sce.all.harmony, group.by = "RNA_snn_res.0.3", features = Lymph, slot = "counts", log = TRUE)

sce.all.harmony@meta.data$cell_type <- sce.all.harmony@meta.data$RNA_snn_res.0.3 %>% as.character %>% as.numeric
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 0] <- "FB"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 1] <- "CM"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 2] <- "CM"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 3] <- "EC"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 4] <- "Myeloid"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 5] <- "Peri"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 6] <- "SMC"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 7] <- "Unknown"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 8] <- "Adi"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 9] <- "Lymphoid"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 10]<- "Neural"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 11]<- "Meso"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 12]<- "FB"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 13]<- "Myeloid"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 14]<- "EC"

p1 <- DimPlot(sce.all.harmony, reduction = "umap",
              group.by = "RNA_snn_res.0.3",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(sce.all.harmony, reduction = "umap",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(sce.all.harmony, reduction = "umap",
              group.by = "orig.ident",
              label = TRUE, pt.size = 0.5) 
p1 + p2 + p3
# ggsave(..., width = 18, height = 6)

# check
if(FALSE){
  tmp <- subset(sce.all.harmony, subset = RNA_snn_res.0.3 == 2)
  DimPlot(tmp, reduction = "umap",
          group.by = "RNA_snn_res.0.3",
          label = TRUE, pt.size = 0.5) 
  
  colors <- c("red3","grey","grey","grey","grey","grey",
              "grey","grey","grey","grey","grey",
              "grey","grey","grey","grey")
  DimPlot(sce.all.harmony, reduction = "umap",
          group.by = "RNA_snn_res.0.3", cols = colors,
          label = TRUE, pt.size = 0.5) 
}

save(sce.all.harmony, file = "sce.all.harmony.annot.Rdata")
# pseudobulk--------------------------------------------------------------------
sce.CM <- sce.all.harmony %>% subset(subset = cell_type == "CM")
sce.EC <- sce.all.harmony %>% subset(subset = cell_type == "EC")
sce.FB <- sce.all.harmony %>% subset(subset = cell_type == "FB")

CM <- sce.CM@assays$RNA$counts %>% as.data.frame %>% 
  t %>% as.data.frame %>% 
  mutate(sample = sce.CM@meta.data$orig.ident) %>% 
  pivot_longer(cols = rownames(sce.CM@assays$RNA$counts), names_to = "symbol", values_to = "expression") %>% 
  group_by(sample, symbol) %>% 
  mutate(sum = sum(expression)) %>% 
  dplyr::select(-expression) %>% 
  filter(!duplicated(paste0(sample,"+",symbol))) %>% 
  pivot_wider(names_from = symbol, values_from = sum) %>% 
  tibble::column_to_rownames("sample")

EC <- sce.EC@assays$RNA$counts %>% as.data.frame %>% 
  t %>% as.data.frame %>% 
  mutate(sample = sce.EC@meta.data$orig.ident) %>% 
  pivot_longer(cols = rownames(sce.EC@assays$RNA$counts), names_to = "symbol", values_to = "expression") %>% 
  group_by(sample, symbol) %>% 
  mutate(sum = sum(expression)) %>% 
  dplyr::select(-expression) %>% 
  filter(!duplicated(paste0(sample,"+",symbol))) %>% 
  pivot_wider(names_from = symbol, values_from = sum) %>% 
  tibble::column_to_rownames("sample")

FB <- sce.FB@assays$RNA$counts %>% as.data.frame %>% 
  t %>% as.data.frame %>% 
  mutate(sample = sce.FB@meta.data$orig.ident) %>% 
  pivot_longer(cols = rownames(sce.FB@assays$RNA$counts), names_to = "symbol", values_to = "expression") %>% 
  group_by(sample, symbol) %>% 
  mutate(sum = sum(expression)) %>% 
  dplyr::select(-expression) %>% 
  filter(!duplicated(paste0(sample,"+",symbol))) %>% 
  pivot_wider(names_from = symbol, values_from = sum) %>% 
  tibble::column_to_rownames("sample")

save(CM, EC, FB, file = "pseudobulk_CM_EC_FB.Rdata")

# Differential Analysis---------------------------------------------------------
# Enrichment Analysis-----------------------------------------------------------
# ORA enrichment----------------------------------------------------------------
# GSEA enrichment---------------------------------------------------------------