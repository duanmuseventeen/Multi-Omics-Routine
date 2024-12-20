# Decoding the transcriptome of calcified atherosclerotic plaque at single-cell resolution
# PMID: 36224302
# 10.1038/s42003-022-04056-7

GSM4837523	Patient 1 AC scRNA-seq
GSM4837524	Patient 1 PA scRNA-seq
GSM4837525	Patient 2 AC scRNA-seq
GSM4837526	Patient 2 PA scRNA-seq
GSM4837527	Patient 3 AC scRNA-seq
GSM4837528	Patient 3 PA scRNA-seq

The cell set samples are indicated by a dash number in the barcodes of the aggregated cell set. The following table shows the dash numbers that correspond to each sample.

Patient 1 PA scRNA-seq barcode -1 suffix
Patient 1 AC scRNA-seq barcode -2 suffix
Patient 2 PA scRNA-seq barcode -3 suffix
Patient 2 AC scRNA-seq barcode -4 suffix
Patient 3 PA scRNA-seq barcode -5 suffix
Patient 3 AC scRNA-seq barcode -6 suffix

# Raw single-cell sequencing data were processed independently per sample using the 
# 10X Genomic Chromium platform77. 10X Genomic’s Cell Ranger with default parameters 
# for read mapping80, sample quality control, unique molecular identifier filtration, 
# normalization, and expression quantification. Further quality filtering is performed 
# by Seurat79 to remove 
# cells with >10% mitochondrial mRNA (total mRNA), 
# in addition to cells with <200 or >4000 genes expressed. 
# Quality control resulted in a reduction of 6145 cells, down to 45,836 cells total
# prior to downsampling. Cell sets were then down-sampled to 17,100 cells to account for 
# imbalances in total cell counts across samples which may adversely influence 
# clustering analysis.

# Doublet analysis and filtering
# In order to identify and remove distressed cells, or rare cell types and artifacts, 
# we first estimated doublet-rates and attempted to filter doublets using the Scrublet 
# package82. However, we found that Scrublet and other doublet detection packages were 
# not effective at identifying doublets in highly heterogenous samples like 
# atherosclerotic plaque. 
# Therefore, as an alternative we identify doublets based on a combination of the 
# inappropriate expression of cell-type marker genes coupled with elevated read counts. 
# Marker genes that should be ubiquitously expressed (>90%) in one cell type (partition) 
# and rarely expressed (<10%) in other cell types (partitions) were used to mark potential 
# doublet cells. Cells with multiple inappropriate marker genes expressed (≥2) are 
# tagged for removal. For example, CD2 expression (cell adhesion molecule specific 
# to T- and NKT-cells) was detected at higher-than-expected levels in VSMCs (3%) and 
# ECs (2.1%), and thus these cells were excluded from analysis. An upward shift in 
# read counts across all doublets relative to the population average was used to 
# validate the doublet identification strategy. In addition, we perform gene co-expression 
# network reconstruction using a custom partial correlations approach (described below), 
# which is sensitive to rare outlying correlation events, and enables us to detect 
# inappropriate co-expression of genes. Inappropriate networks reflective of unexpected 
# cell-types were identified during the doublet detection process to validate our 
# filtration strategy. Remaining cells were re-clustered producing six major partitions; 
# macrophages, ECs, VSMCs, NKT cells, T- and B-lymphocytes (Fig. 1c, d). Patient-specific 
# cells are shown in Fig. S3a.

rm(list = ls());gc()
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree) # use for determine the optimal resolution
library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
library(Matrix)
library(scDblFinder)
library(DoubletFinder)
library(sctransform)
# Load data---------------------------------------------------------------------
# h5
if(TRUE){
  # hub <- c(
  #   "GSM4837523_02dat20190515tisCARconDIS_moleculeinfo.h5",               
  #   "GSM4837524_01dat20190515tisCARconHEA_moleculeinfo.h5",               
  #   "GSM4837525_02dat20190620tisCARconDIS_moleculeinfo.h5",               
  #   "GSM4837526_01dat20190620tisCARconHEA_moleculeinfo.h5",               
  #   "GSM4837527_02dat20190717tisCARconDIS_moleculeinfo.h5",               
  #   "GSM4837528_01dat20190717tisCARconHEA_moleculeinfo.h5"    
  # )
  # 
  # scobjList <- lapply(hub, function(x){
  #   sce <- x %>% 
  #     Read10X_h5 %>% 
  #     CreateSeuratObject(project = x %>% str_extract("^GSM[0-9]*"))
  #   return(sce)
  # })

  # 错误于x$exists(name): HDF5-API Errors:
  #   error #000: ../../src/H5L.c in H5Lexists(): line 943: unable to get link info
  # class: HDF5
  # major: Links
  # minor: Can't get value
  # 
  #   error #001: ../../src/H5VLcallback.c in H5VL_link_specific(): line 5173: unable to execute link specific callback
  #       class: HDF5
  #       major: Virtual Object Layer
  #       minor: Can't operate on object
  # 
  # error #002: ../../src/H5VLcallback.c in H5VL__link_specific(): line 5136: unable to execute link specific callback
  # class: HDF5
  # major: Virtual Object Layer
  # minor: Can't operate on object
  # 
  #   error #003: ../../src/H5VLnative_link.c in H5VL__native_link_specific(): line 329: unable to specific link info
  #       class: HDF5
  #       major: Links
  #       minor: Object not found
  # 
  #   error #004: ../../src/H5L.c in H5L__exists(): line 3082: path doesn't exist
  # class: HDF5
  # major: Links
  # minor: Object already exists
  # 
  # error #005: ../../src/H5Gtraverse.c 
  
  # scobj <- merge(x = scobjList[[1]], y = scobjList[-1], add.cell.ids = dir())
}
# mtx
if(TRUE){
  hub <- c(
    "GSM4837523",                                                         
    "GSM4837524",                                                         
    "GSM4837525",                                                         
    "GSM4837526",                                                         
    "GSM4837527",                                                         
    "GSM4837528"  
  )
  names(hub) <- hub
  
  scobj <- hub %>% 
    Read10X %>% 
    CreateSeuratObject(
      project = "GSE159677",
      min.cells = 0,
      min.features = 0,
      assay = "RNA")                                              
}

# head(scobj@meta.data)
# orig.ident nCount_RNA nFeature_RNA
# GSM4837523_AAACCCAAGATTAGAC-1 GSM4837523      11667         2905
# GSM4837523_AAACCCAAGCATGTTC-1 GSM4837523       5385         1494
# GSM4837523_AAACCCAAGCCTGTCG-1 GSM4837523       1471          575
# GSM4837523_AAACCCAAGGGTTTCT-1 GSM4837523       8957         2515
# GSM4837523_AAACCCAAGTGCACCC-1 GSM4837523      11681         2586
# GSM4837523_AAACCCACAAGTTTGC-1 GSM4837523       2495          794

# table (scobj@meta.data$orig.ident)
# GSM4837523 GSM4837524 GSM4837525 GSM4837526 GSM4837527 GSM4837528 
# 11015       3716      15960       5523      12388       3379 
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

# VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
# FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# data.frame(
#   nFeature = scobj@meta.data$nFeature_RNA,
#   group = scobj@meta.data$orig.ident
# ) %>% 
#   ggplot(aes(x = nFeature, color = group)) +
#   geom_density() + 
#   scale_x_continuous(breaks = c(200,300,400,500,750,1000,3000,4000,5000))

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 750 &
    nFeature_RNA < 4000 & 
    nCount_RNA < 20000 &
    percent.mt < 10 &
    percent.rp < 40 &
    percent.hb < 0.1)  

# VlnPlot(scobj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)

# down sample
# https://cloud.tencent.com/developer/article/1825648
table (scobj@meta.data$orig.ident)
# GSM4837523 GSM4837524 GSM4837525 GSM4837526 GSM4837527 GSM4837528 
# 8754       2937      12231       4274       9029       2536 

scobj.ds <- subset(scobj, downsample = 2500)

table (scobj.ds@meta.data$orig.ident)
# GSM4837523 GSM4837524 GSM4837525 GSM4837526 GSM4837527 GSM4837528 
# 2500       2500       2500       2500       2500       2500 

# scobj.ds <- scobj.ds %>% 
#   NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
#   ScaleData %>% 
#   RunPCA(npcs = 50)

scobj.ds <- scobj.ds %>% 
  SCTransform %>% # vst.flavor = 'v2', verbose = FALSE
  RunPCA(npcs = 50)
# Integration-------------------------------------------------------------------
scobj.ds.harmony <- RunHarmony(
  scobj.ds,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  reduction.save = "harmony")
# Cluster-----------------------------------------------------------------------
scobj.ds.harmony <- FindNeighbors(scobj.ds.harmony, reduction = "harmony", dims = 1:30)
scobj.ds.harmony <- FindClusters(scobj.ds.harmony, resolution = seq(0.1, 1, 0.1))

clustree(scobj.ds.harmony, prefix = "SCT_snn_res.")
# Visualization-----------------------------------------------------------------
#UMAP
scobj.ds.harmony <- RunUMAP(scobj.ds.harmony, reduction = "harmony", min_dist = 0.3, dims = 1:30)
#T-SNE
scobj.ds.harmony <- RunTSNE(scobj.ds.harmony, reduction = "harmony", dims = 1:30)

# Here we examine the single-cell transcriptome of entire calcified atherosclerotic 
# core (AC) plaques and patient-matched proximal adjacent (PA) portions of carotid 
# artery tissue from patients undergoing carotid endarterectomy. 

meta <- scobj.ds.harmony@meta.data
meta <- meta %>% 
  mutate(disease = case_when(orig.ident %in% c("GSM4837523","GSM4837525","GSM4837527") ~ "AC",
                             orig.ident %in% c("GSM4837524","GSM4837526","GSM4837528") ~ "PA"))
scobj.ds.harmony@meta.data <- meta

DimPlot(scobj.ds.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.ds.harmony, group.by = "disease", reduction = "umap", label = T)
DimPlot(scobj.ds.harmony, group.by = "SCT_snn_res.0.4", reduction = "umap", label = T)
DimPlot(scobj.ds.harmony, group.by = "SCT_snn_res.0.4", reduction = "tsne", label = T)
# Cell cycle--------------------------------------------------------------------
exp.mat <- read.table(file = "/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.ds.harmony <- CellCycleScoring(scobj.ds.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(scobj.ds.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

scobj.ds.harmony <- RunPCA(scobj.ds.harmony, features = c(s.genes, g2m.genes))
DimPlot(scobj.ds.harmony, reduction = "pca")
DimPlot(scobj.ds.harmony, reduction = "umap")
# Filter doublets---------------------------------------------------------------
if(TRUE){
  # scDblFinder
  scobj.ds.harmony.sce <- as.SingleCellExperiment(scobj.ds.harmony, assay = "SCT")
  scobj.ds.harmony.sce <- scDblFinder(scobj.ds.harmony.sce, dbr.sd=1)
  
  scobj.ds.harmony <- as.Seurat(scobj.ds.harmony.sce)
  scobj.h.sc <- subset(scobj.ds.harmony, scDblFinder.class == "singlet")
}
# Annotation--------------------------------------------------------------------
scobj.h.sc <- RegroupIdents(scobj.h.sc, "SCT_snn_res.0.4")

# 36224302
if(TRUE){
  Macrophage <- c("AIF1","CD14","CD68") # 4 6 7 ?14 ?17
  Endothelial <-c("VWF","PECAM1","ECSCR") # 0 10
  VSMC <- c("CALD1","MYL9","TAGLN") # 2 5 8 ?10 13 
  NKT <- c("NKG7","XCL1","CTSW","CD69") # 1 3 9 12 ?15 16
  T_cell <-c("CD2","TRAC","CD69","CD3E","CD4","CD8A") # 1 3 9 12 16
  B_cell <-c("CD69","CD79A","MS4A1","IGKC","CD19") # 11
}
# 36172868
if(TRUE){
  SMC <- c('MYH11', 'ACTA2','CNN1', 'SMTN', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') # 2 8 13
  Peri <- c('RGS5', 'ABCC9', 'KCNJ8') # ?2 ?5 13
  FB <- c('LAMA2', 'TSHZ2') # 5
  EC <- c('PECAM1') 
  RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ')
  Meso <- c('MSLN', 'WT1', 'BNC1') 
  Adi <- c('GPAM', 'FASN', 'LEP')
  Neural <- c('PLP1', 'NRXN1', 'NRXN3') 
  Lymph <- c('CD3E', 'IL7R', 'CD40LG') 
  M <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74', 'CD14', 'C1QA') 
}
DC <- c("FCER1A", "CD1E", "CD1C", "HLA-DMA", "HLA-DMB") # 

FindMarkers(scobj.h.sc, ident.1 = 17, group.by = "SCT_snn_res.0.4", min.pct = 0.25)

VlnPlot(scobj.h.sc, features = DC, group.by = "SCT_snn_res.0.4", layer = "counts", log = TRUE)

scobj.h.sc@meta.data$cell_type <- "Unknown"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.4 %in% c(0,10)] <- "EC"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.4 %in% c(1,3,9,12,15,16)] <- "Lymphoid"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.4 %in% c(2,5,8,13)] <- "VSMC"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.4 %in% c(11)] <- "B"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.4 %in% c(4,6,7,14,17)] <- "Myeloid"

p1 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "SCT_snn_res.0.4",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "orig.ident",
              label = TRUE, pt.size = 0.5) 
p4 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "disease",
              label = TRUE, pt.size = 0.5) 
p_out <- p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2)

p5 <- RidgePlot(object = scobj.h.sc, group.by = "cell_type", features = "NEU1")
p6 <- VlnPlot(object = scobj.h.sc, group.by = "cell_type", features = "NEU1")

p7 <- RidgePlot(object = scobj.h.sc, group.by = "cell_type", features = "COL1A1")
p8 <- VlnPlot(object = scobj.h.sc, group.by = "cell_type", features = "COL1A1")

save(scobj.h.sc, file = "GSE159677-scobj.h.sc.annot.Rdata")
ggsave(p_out, filename = "GSE1159677.pdf", width=10, height=10, units="in")
ggsave(p5, filename = "GSE115469-Ridge-NEU1.pdf", width=5, height=5, units="in")
ggsave(p6, filename = "GSE115469-Vln-NEU1.pdf", width=5, height=5, units="in")
ggsave(p7, filename = "GSE115469-Ridge-COL1A1.pdf", width=5, height=5, units="in")
ggsave(p8, filename = "GSE115469-Vln-COL1A1.pdf", width=5, height=5, units="in")
