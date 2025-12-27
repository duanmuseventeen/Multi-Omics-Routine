# Resolving the intertwining of inflammation and fibrosis in human heart failure at single‑cell level----
# GSE145154

# Resolving the intertwining of inflammation and fibrosis in human heart failure at single-cell level
# Inflammation and fibrosis are intertwined mechanisms fundamentally involved in 
# heart failure. Detailed deciphering gene expression perturbations and cell–cell 
# interactions of leukocytes and non-myocytes is required to understand cell-type-specific 
# pathology in the failing human myocardium. To this end, we performed single-cell 
# RNA sequencing and single T cell receptor sequencing of 200,615 cells in both human 
# dilated cardiomyopathy (DCM) and ischemic cardiomyopathy (ICM) hearts. We sampled both 
# lesion and mild-lesion tissues from each heart to sequentially capture cellular and 
# molecular alterations to different extents of cardiac fibrosis. By which, left (lesion) 
# and right ventricle (mild-lesion) for DCM hearts were harvest while infarcted (lesion) 
# and non-infarcted area (mild-lesion) were dissected from ICM hearts. A novel 
# transcription factor AEBP1 was identified as a crucial cardiac fibrosis regulator 
# in ACTA2+ myofibroblasts. Within fibrotic myocardium, an infiltration of a 
# considerable number of leukocytes was witnessed, especially cytotoxic and exhausted 
# CD8+ T cells and pro-inflammatory CD4+ T cells. Furthermore, a subset of tissue-resident 
# macrophage, CXCL8hiCCR2+HLA-DRhi macrophage was particularly identified in severely 
# fibrotic area, which interacted with activated endothelial cell via DARC, that 
# potentially facilitate leukocyte recruitment and infiltration in human heart failure.

# First, we split the combined object into a list using “SplitObject” function,
# with each donor as an element. Prior to finding anchors, we performed
# log-normalization and identified 2000 variable features individually for each
# with the “vst” method.
# 
# Seurat::SplitObject
# 
# scobj <- scobj %>% 
#   NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# 
# Next, we identified anchors using the “FindIntegrationAnchors” function with 
# default parameters. Here, we specified the normal datasets as a “reference”, with
# the remainder designated as “query” datasets. This approach substantially reduced 
# computational time, particularly a large number of datasets. 
# We then integrated the batches using the anchors with the “IntegrateData” function with
# the “dims” parameter set to 30, which returns a Seurat object with a batch-corrected 
# expression matrix for all cells. 
# 
# ?Seurat::FindIntegrationAnchors
# ?Seurat::IntegrateData
# 
# We then used this new integrated matrix for downstream analysis and visualization.
# We scaled the integrated data, setting the parameter “vars.to.regress” to 
# “percent.mito” and “nCount_RNA”. Principal component analysis (PCA) was performed 
# using the “RunPCA” function.
# 
# %>% 
#   ScaleData %>% 
#   RunPCA(npcs = 50)


rm(list = ls())

# Load pkgs---------------------------------------------------------------------
require(dplyr)
require(Seurat)
require(harmony)
require(ggplot2)
require(clustree)
require(ROGUE)
require(stringr)
# Prepare-----------------------------------------------------------------------
for(i in dir()){
  setwd(paste0("GSE145154/scRNAseq/", i))
  file.rename(dir(), str_remove(dir(), "^GSM[0-9A-Z_-]*"))
}
# Load data---------------------------------------------------------------------
hub <- c("DCM_2_LVN",  
         "DCM_2_RVP",
         "DCM_3_RVN",
         "ICM_1_MIN",   
         "ICM_1_NMIP",  
         "ICM_2_RVN",  
         "ICM_3_LVP",  
         "N_1_LVN",  
         "N_1_RVP",
         "DCM_2_LVP",  
         "DCM_3_LVN",  
         "DCM_3_RVP",              
         "ICM_1_MIP",   
         "ICM_2_LVN",   
         "ICM_2_RVP",  
         "ICM_3_RVN",
         "N_1_LVP",
         "DCM_2_RVN",  
         "DCM_3_LVP",
         "ICM_1_NMIN",
         "ICM_2_LVP",
         "ICM_3_LVN",
         "ICM_3_RVP",
         "N_1_RVN")
names(hub) <- hub

scobj <- hub %>% 
  Read10X %>%  
  CreateSeuratObject(project = "GSE145154")
# add meta.data-----------------------------------------------------------------
scobj@meta.data$disease <- scobj@meta.data$orig.ident 
scobj@meta.data$orig.ident <- str_remove_all(rownames(scobj@meta.data), "_[A-Z]*_[A-Z0-9-]*$")
scobj@meta.data$tissue <- rownames(scobj@meta.data) %>% str_remove_all("_[A-Z0-1-]*$") %>% str_remove_all("^[A-Z]*_[0-9]_")
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")

pdf("QC.pdf")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Cells were discarded according to the following
# criteria: (1) cells that had fewer than 500 genes (UMI > 0);
# (2) cells that had fewer than 800 UMI or over 8000 UMI;
# and (3) cells that had more than 10% mitochondrial UMI counts.
scobj <- subset(scobj, subset = nCount_RNA >= 500 & nFeature_RNA > 800 & nFeature_RNA < 8000 & percent.mt < 10) 
  
pdf("QC_2.pdf")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)
FeatureScatter(scobj, 
               slot = "counts",
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA")
dev.off()

scobj <- scobj %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>%
  RunPCA(npcs = 50)

pdf("PCA_QC.pdf")
DimHeatmap(scobj, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()
# Integrate data----------------------------------------------------------------
scobj <- RunHarmony(scobj,
                    group.by.vars = "orig.ident",
                    reduction.use = "pca",
                    reduction.save = "harmony")

scobj <- FindNeighbors(scobj, reduction = "harmony", dims = 1:30)
scobj <- FindClusters(scobj, resolution = seq(0.1, 1, 0.1))

# select optimal res # 0.2 0.3 0.4 0.5
# 1 clustree
pdf("hsample_clustree.pdf")
clustree(scobj, prefix = "RNA_snn_res.")
dev.off()
# # 2 marker gene AUC
# scobj %>% 
#   RegroupIdents(metadata = "RNA_snn_res.0.2") %>% 
#   FindAllMarkers(test.use = 'roc') %>% 
#   filter(myAuc > 0.6) %>% 
#   count(cluster, name = "number")
# 
# # 3 ROGUE
# pdf("ROGUE.pdf")
# 
# scobj_sub <- scobj %>% subset(downsample = 1000)
# expr <- GetAssayData(scobj_sub, layer = "counts") %>% as.matrix
# rogue(expr, labels = scobj_sub@meta.data$RNA_snn_res.0.2,
#       samples = scobj_sub@meta.data$orig.ident, platform='UMI', span =0.6) %>%
#   rogue.boxplot()
# 
# rogue(expr, labels = scobj_sub@meta.data$RNA_snn_res.0.4,
#       samples = scobj_sub@meta.data$orig.ident, platform='UMI', span=0.6) %>%
#   rogue.boxplot()
# 
# dev.off()
# Visualization-----------------------------------------------------------------
#UMAP
scobj <- RunUMAP(scobj, reduction = "harmony", min_dist = 0.3, dims = 1:30)
#T-SNE
scobj <- RunTSNE(scobj, reduction = "harmony", dims = 1:30)

pdf("select resolution.pdf")
scobj <- RegroupIdents(scobj, metadata = "RNA_snn_res.0.2")
DimPlot(scobj, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj, group.by = "RNA_snn_res.0.2", reduction = "umap", label = T)
DimPlot(scobj, group.by = "orig.ident", reduction = "tsne", label = T)
DimPlot(scobj, group.by = "RNA_snn_res.0.2", reduction = "tsne", label = T)

scobj <- RegroupIdents(scobj, metadata = "RNA_snn_res.0.5")
DimPlot(scobj, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj, group.by = "RNA_snn_res.0.5", reduction = "umap", label = T)
DimPlot(scobj, group.by = "orig.ident", reduction = "tsne", label = T)
DimPlot(scobj, group.by = "RNA_snn_res.0.5", reduction = "tsne", label = T)
dev.off()
# Annotation------------------------------------------------------------------------------
scobj <- RegroupIdents(scobj, metadata = "RNA_snn_res.0.5")
ACM <- c('TNNT2', 'TTN', 'MYH6', 'MYH7') # 12
VCM <- c('MYH7', 'MYL2', 'FHL2') # 12 
EC <- c('VWF', 'CDH5', 'THBD', 'TEK', 'ENG', 'PECAM1') # 5 10 11 14 16 18
FB <- c('FN1', 'VIM', 'COL1A1', 'COL1A2') # 2 ?7 ?22
M <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74') # 1 4 8 13 ?19 ?20
Myeloid	<- c('CD163', 'CD14', 'C1QA') # 1 4 8 13 ?19
P <- c('HCN1', 'HCN4', 'CACNA1D') # pacemaker cells (P cells)
Peri <- c('RGS5', 'ABCC9', 'KCNJ8') # 7 14
SMC <- c('MYH11', 'ACTA2','CNN1', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') # ?7 9 ?14
Meso <- c('MSLN', 'WT1', 'BNC1') # ?18 
Neural <- c('PLP1', 'NRXN1', 'NRXN3') # 22
Adi <- c('GPAM', 'FASN', 'LEP') # 
Lymph <- c('CD3E', 'IL7R', 'CD40LG') # 0 ?3 ?6 ?13 ?17 
Mast <- c('KIT', 'CPA3') # 21
B <- c("CD19", "CD20") # 15
RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ')

if(TRUE){
  VlnPlot(scobj, features = ACM, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = VCM, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = EC, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = FB, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = M, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = Myeloid, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = P, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = Peri, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = SMC, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = Meso, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = Neural, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = Adi, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = Lymph, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = Mast, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = B, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
  VlnPlot(scobj, features = RBC, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)
}

scobj@meta.data$cell_type <- "Unknown"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.5 %in% c(12)] <- "CM"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.5 %in% c(5, 10, 11, 14, 16, 18)] <- "EC"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.5 %in% c(2)] <- "FB"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.5 %in% c(1, 4, 8, 13, 15, 19, 20)] <- "Myeloid"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.5 %in% c(7)] <- "Peri"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.5 %in% c(9)] <- "SMC"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.5 %in% c(22)] <- "Neural"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.5 %in% c(21)] <- "Mast"
scobj@meta.data$cell_type[scobj@meta.data$RNA_snn_res.0.5 %in% c(0, 3, 6, 17)] <- "Lymphoid"

p1 <- DimPlot(scobj, reduction = "umap",
              group.by = "RNA_snn_res.0.5",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj, reduction = "umap",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj, reduction = "umap",
              group.by = "orig.ident",
              label = TRUE, pt.size = 0.5) 
p1 + p2 + p3

save(scobj.annot, file = "GSE145154/GSE145154-scobj.annot.Rdata")
