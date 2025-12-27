# Single cell and lineage tracing studies reveal the impact of CD34+ cells on myocardial fibrosis during heart failure
# GSE222144
# https://explore.data.humancellatlas.org/projects/86fe0a0c-88b3-4a3e-94a1-6f9feadc401e
# scRNA-seq for total non-cardiomyocytes from human heart

# Reference
# https://blog.csdn.net/lijianpeng0302/article/details/130892848
# https://zhuanlan.zhihu.com/p/562315878
# Automatic cell type identification methods for single-cell RNA sequencing

rm(list = ls())
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree) # use for determine the optimal resolution
library(ROGUE) # use for determine the optimal resolution
library(harmony)
# Load data---------------------------------------------------------------------
import <- function(dir, project){
  res <- Read10X(data.dir = dir) |>
    CreateSeuratObject(project = project, min.cells = 3, min.features = 200)
  return(res)
}

Control <- import(dir = "GSE222144/Control/", project = "Ctrl")
HF <- import(dir = "GSE222144/HF/", project = "HF")

scobj <- merge(x = Control, y = HF, project = "GSE222144", add.cell.ids = c("Ctrl","HF"))
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

if(TRUE){
  data.frame(
    nFeature = c(Control@meta.data$nFeature_RNA,
                 HF@meta.data$nFeature_RNA),
    group = c(rep("Control", length(Control@meta.data$nFeature_RNA)), 
              rep("HF", length(HF@meta.data$nFeature_RNA)))
  ) %>% 
    ggplot(aes(x = nFeature, color = group)) +
    geom_density() + 
    scale_x_continuous(breaks = c(200,500,1000,3000,4000,5000))
}

scobj <- subset(scobj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10 & percent.rp < 20)  
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)

scobj <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50)
# Integration-------------------------------------------------------------------
# Perform integration (harmony)
scobj <- IntegrateLayers(
    object = scobj,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = 'harmony',
    verbose = FALSE)

# scobj <- IntegrateLayers(
#   object = scobj, 
#   method = CCAIntegration, 
#   orig.reduction = "pca", 
#   new.reduction = 'integrated.cca',
#   verbose = FALSE)

# re-join layers after integration
scobj[["RNA"]] <- JoinLayers(scobj[["RNA"]])
# Cluster-----------------------------------------------------------------------
scobj <- FindNeighbors(scobj, reduction = "harmony", dims = 1:30)
scobj <- FindClusters(scobj, resolution = seq(0.1, 1, 0.1))
clustree(scobj, prefix = "RNA_snn_res.")
# Visualization-----------------------------------------------------------------
#UMAP
scobj <- RunUMAP(scobj, reduction = "harmony", min_dist = 0.3, dims = 1:30)
# head(sce.all@reductions$umap@cell.embeddings) 
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(scobj, group.by = "RNA_snn_res.0.4", reduction = "umap", label = T)
#T-SNE
scobj <- RunTSNE(scobj, reduction = "harmony", dims = 1:30)
# head(sce.all@reductions$tsne@cell.embeddings)
DimPlot(scobj, group.by = "RNA_snn_res.0.4", reduction = "tsne", label = T) 

# Annotation--------------------------------------------------------------------
scobj <- RegroupIdents(scobj, metadata = "RNA_snn_res.0.4")
ACM <- c('TNNT2', 'TTN', 'MYH6', 'MYH7') # ?6 8
VCM <- c('MYH7', 'MYL2', 'FHL2') # 8
EC <- c('VWF', 'CDH5', 'THBD', 'TEK', 'ENG', 'PECAM1') # 1 ?4 11
FB <- c('FN1', 'VIM', 'COL1A1', 'COL1A2') # 0 9
M <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74') # 5 6
Myeloid	<- c('CD163', 'CD14', 'C1QA') # 6
P <- c('HCN1', 'HCN4', 'CACNA1D') # pacemaker cells (P cells) ?3
Peri <- c('RGS5', 'ABCC9', 'KCNJ8') # 3
SMC <- c('MYH11', 'ACTA2','CNN1', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') # 3
Meso <- c('MSLN', 'WT1', 'BNC1') # 
Neural <- c('PLP1', 'NRXN1', 'NRXN3') # 10
Adi <- c('GPAM', 'FASN', 'LEP') # 
Lymph <- c('CD3E', 'IL7R', 'CD40LG') # 7
Mast <- c('KIT', 'CPA3') # 
RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ') # ?12

VlnPlot(scobj, features = EC, group.by = "RNA_snn_res.0.4", layer = "counts", log = TRUE)

scobj@meta.data$cell_type <- 0
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 0] <- "FB"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 1] <- "EC"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 2] <- "Unknown"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 3] <- "SMC"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 4] <- "EC"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 5] <- "Myeloid" # need to subtype
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 6] <- "Myeloid" # need to subtype
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 7] <- "Lymphoid" # need to subtype
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 8] <- "CM"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 9] <- "FB"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 10]<- "Neural"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 11]<- "EC"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.4 == 12]<- "HBC"

p1 <- DimPlot(scobj, reduction = "umap",
              group.by = "RNA_snn_res.0.4",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj, reduction = "umap",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj, reduction = "umap",
              group.by = "orig.ident",
              label = TRUE, pt.size = 0.5) 
p1 + p2 + p3
# ggsave(..., width = 18, height = 6)

markers.2 <- FindMarkers(scobj, group.by = "RNA_snn_res.0.4", ident.1 = "2", min.pct = 0.25)
markers.4 <- FindMarkers(scobj, group.by = "RNA_snn_res.0.4", ident.1 = "4", min.pct = 0.25)

markers.2 # from 0: FB
markers.4 # from 1: EC

# save(scobj, file = "GSE222144/GSE222144-scobj.annot.Rdata")
load("GSE222144/GSE222144-scobj.Rdata")
# Auto--------------------------------------------------------------------------
## √ SingleR-----------------------------------------------------------------------
library(SingleR)
library(celldex)

hpca.se <- celldex::HumanPrimaryCellAtlasData()
hpca.se
hpca.se@metadata

# bpe.se<-celldex::BlueprintEncodeData() 
# bpe.se@metadata

pred.scobj <- SingleR(
  test = GetAssayData(scobj, layer = "data"), 
  ref = hpca.se, 
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  labels = hpca.se$label.main)

all(pred.scobj@rownames == colnames(scobj@assays$RNA$data))
[1] TRUE
table(pred.scobj$labels, scobj@meta.data$cell_type)

scobj@meta.data$singleR <- pred.scobj$labels
p4 <- DimPlot(scobj, reduction = "umap",
        group.by = "singleR",
        label = TRUE, pt.size = 0.5) 
p2 + p4
## * scmap-------------------------------------------------------------------------
library(SingleCellExperiment)
library(scmap)
## CaSTLe-------------------------------------------------------------------------
## * CHETAH-------------------------------------------------------------------------
## * clustifyr-------------------------------------------------------------------------
## Garnett-------------------------------------------------------------------------
## scClassifR-------------------------------------------------------------------------
## SciBet-------------------------------------------------------------------------
## scID-------------------------------------------------------------------------
## scLearn-------------------------------------------------------------------------
## scmap-cluster-------------------------------------------------------------------------
## scmapcell-------------------------------------------------------------------------
## * [笔记本带不动]scPred------------------------------------------------------------------------
# https://github.com/powellgenomicslab/scPred
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1862-5

devtools::install_github("powellgenomicslab/scPred")
library(scPred)                     
         
sce.all <- readRDS("Defining cardiac functional recovery in end-stage heart failure at single-cell resolution/GSE226314_global.rds")
sce.all <- sce.all %>% subset(subset = condition %in% c("Donor","Rpre","NRpre"))
reference <- getFeatureSpace(sce.all, "cell.type")
reference <- trainModel(reference)

scobj <- scPredict(scobj, reference)
DimPlot(scobj, group.by = "scpred_prediction", reduction = "scpred")

scobj <- RunUMAP(scobj, reduction = "scpred", dims = 1:15)
DimPlot(scobj, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
## * [运行出错]Seurat------------------------------------------------------------------------
library(SeuratData)

# select two technologies for the query datasets
AD <- AvailableData()
AD %>% filter(species == "human" & system == "heart")

# InstallData("heartref", force.reinstall = TRUE)

# install.packages(
#   "GSE222144/heartref.SeuratData_1.0.0.tar.gz",
#   repos = NULL,
#   type = "source")
# 
# library(heartref.SeuratData)
# LoadData("heartref")

heartref <- readRDS("GSE222144/heartref.SeuratData_1.0.0/ref.Rds")
heartref <- NormalizeData(heartref, normalization.method = "LogNormalize", scale.factor = 10000)
错误于validObject(object = value): 
  类别为"SCTAssay"的对象无效: 'meta.features' must have the same number of rows as 'data'

# heartref <- FindTransferAnchors(
#   reference = heartref, 
#   query = heartref, 
#   dims = 1:30,
#   reference.reduction = "pca")
# predictions <- TransferData(anchorset = heartref.anchors, refdata = heartref.ref$celltype, dims = 1:30)
# heartref.query <- AddMetaData(heartref.query, metadata = predictions)
# 
# heartref.query$prediction.match <- heartref.query$predicted.id == heartref.query$celltype
# table(heartref.query$prediction.match)
# 
# table(heartref.query$predicted.id)
# VlnPlot(heartref.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")
# 
# heartref.ref <- RunUMAP(heartref.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
# heartref.query <- MapQuery(anchorset = heartref.anchors, reference = heartref.ref, query = heartref.query,
#                            refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
# 
# p1 <- DimPlot(heartref.ref, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
#               repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
# p2 <- DimPlot(heartref.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
#               label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
# p1 + p2

## SingleCellNet-------------------------------------------------------------------------
## *CellAssign--------------------------------------------------------------------
# https://zhuanlan.zhihu.com/p/88924349
## * scCATCH-----------------------------------------------------------------------
## * SCINA-------------------------------------------------------------------------
## scTyper-------------------------------------------------------------------------
## SVM---------------------------------------------------------------------------
## ACTINN-------------------------------------------------------------------------