setwd("GSE193337_NAT vs PCa/")

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
library(DESeq2)
require(stringr)
require(Matrix)
library(copykat)

myqc4seurat <- function(seurat.obj,
                        xintercept1 = c(200,300,400,500,1000,5000,6000,7000,8000),
                        xintercept2 = c(200,500,1000,5000,10000,15000)){
  p1 <- VlnPlot(seurat.obj, features = c("nFeature_RNA"), group.by = "orig.ident", ncol = 1) + scale_y_continuous(breaks = c(200,500, 1000,2000,4000,6000,8000,10000))
  p2 <- VlnPlot(seurat.obj, features = c("nCount_RNA"), group.by = "orig.ident", ncol = 1)
  p3 <- VlnPlot(seurat.obj, features = c("percent.mt"), group.by = "orig.ident", ncol = 1) + scale_y_continuous(breaks = c(10,15,20))
  p4 <- VlnPlot(seurat.obj, features = c("percent.rp"), group.by = "orig.ident",  ncol = 1) + scale_y_continuous(breaks = c(10,20,30,40))
  p5 <- VlnPlot(seurat.obj, features = c("percent.hb"), group.by = "orig.ident", ncol = 1)
  p6 <- FeatureScatter(seurat.obj, group.by = "orig.ident", feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p7 <- data.frame(
    nFeature = seurat.obj$nFeature_RNA,
    group = seurat.obj$orig.ident
  ) %>%
    ggplot(aes(x = nFeature, color = group)) +
    geom_density() +
    geom_vline(xintercept = xintercept1, color = "gray50", linetype = 2) +
    theme_classic()
  
  p8 <-  data.frame(
    nCount = seurat.obj$nCount_RNA,
    group = seurat.obj$orig.ident
  ) %>%
    ggplot(aes(x = nCount, color = group)) +
    geom_density() +
    geom_vline(xintercept = xintercept2, color = "gray50", linetype = 2) +
    theme_classic()
  
  list(p1, p2, p3, p4, p5, p6, p7, p8)
}
# load data---------------------------------------------------------------------
set.seed(1011)

setwd("rawdat/")
hub <- dir()
names(hub) <- hub

scobj <- hub %>%
  Read10X %>%
  CreateSeuratObject(
    min.cells = 0,
    min.features = 0,
    project = "GSE193337",
    assay = "RNA")

scobj
# An object of class Seurat 
# 33538 features across 39516 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
# 1 layer present: counts

setwd("..")
metadata   <- scobj@meta.data
meta       <- readxl::read_excel("GSE193337-GPL20301_series_matrix.xlsx")
meta_merge <- metadata %>% 
  mutate(geo_accession = orig.ident) %>% 
  left_join(meta, by = "geo_accession") %>% 
  as.data.frame
rownames(meta_merge) <- rownames(metadata)
meta_merge -> scobj@meta.data

table(scobj@meta.data$orig.ident)
# GSM5793824 GSM5793825 GSM5793826 GSM5793827 GSM5793828 GSM5793829 GSM5793831 GSM5793832 
# 6674       2644       1875       6557       2574       5922       5358       7912
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

pdf("run1_QC前.pdf")
myqc4seurat(seurat.obj = scobj)
dev.off()

scobj.qc <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    nCount_RNA > 500 &
    nCount_RNA < 15000 &
    percent.mt < 20 &
    percent.rp < 25 &
    percent.hb < 1)  

pdf("run1_QC后.pdf")
myqc4seurat(seurat.obj = scobj.qc)
dev.off()

table(scobj.qc$orig.ident)
# GSM5793824 GSM5793825 GSM5793826 GSM5793827 GSM5793828 GSM5793829 GSM5793831 GSM5793832 
# 4871        971        690       2784       1373       2075       1516       2429
# blacklist---------------------------------------------------------------------
scobj = scobj.qc
scobj@assays$RNA$raw_count = scobj@assays$RNA$counts

blacklist <- readxl::read_excel("blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(scobj))
# [1] 243

scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj)
# [1] 33296 16709
# protein coding----------------------------------------------------------------
gtf_data <- rtracklayer::import("refdata-gex-GRCh38-2024-A_genes.gtf") %>% as.data.frame
keep     <- gtf_data$gene_name[gtf_data$gene_type != "lncRNA"] %>% unique
scobj    <- subset(scobj , features = keep)

dim(scobj)
# [1] 19041 16709
# save -------------------------------------------------------------------------
scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "run1_scobj(GSE193337).Rdata")
# Nomalization & harmony & cluster----------------------------------------------
scobj.harmony <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>%
  JoinLayers(assay = "RNA")

scobj.harmony.20 <- scobj.harmony %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.25 <- scobj.harmony %>%
  FindNeighbors(reduction = "pca", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.30 <- scobj.harmony %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

pdf("run1_nf2000_h clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony, reduction = "pca", ndims = 50)
clustree(scobj.harmony.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.30, prefix = "RNA_snn_res.")
dev.off()

scobj.harmony <- scobj.harmony.30 %>% 
  RunUMAP(reduction = "pca", dims = 1:30)
# Visualization-----------------------------------------------------------------
pdf("run1_nf2000_h_pc20.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "sample", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "group", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.8", reduction = "umap", label = T)
FeaturePlot(scobj.harmony, features = "percent.mt")
dev.off()
# Cell Cycle--------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.harmony <- CellCycleScoring(scobj.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scobj.harmony <- RunPCA(scobj.harmony, features = c(s.genes, g2m.genes), reduction.name = "cellcycle")

pdf("run1_cell cycle.pdf")
RidgePlot(scobj.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
DimPlot(scobj.harmony, reduction = "cellcycle")
DimPlot(scobj.harmony, reduction = "umap")
dev.off()
# doublets----------------------------------------------------------------------
save(scobj.harmony, file = "run1_scobj.harmony(GSE193337).Rdata")

scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list

nPC = 30
for (i in 1:length(scobj.harmony.split)) {
  # pK Identification (ground-truth) -------------------------------------------
  sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:nPC)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
  ## Homotypic Doublet Proportion Estimate -------------------------------------
  homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]@meta.data$RNA_snn_res.0.8)  
  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi <- round(0.05 * nrow(scobj.harmony.split[[i]]@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies ----------------
  # scobj.harmony.split[[i]] <- doubletFinder(scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/228
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/225#issuecomment-2786505997
  scobj.harmony.split[[i]] <- doubletFinder(
    scobj.harmony.split[[i]], PCs = 1:nPC, pN = 0.25, pK = pK, 
    nExp = nExp_poi.adj, sct = FALSE)
}

save(scobj.harmony.split, file = "run1_scobj.harmony.split(GSE193337).Rdata")

Singlet <- c()
for (i in 1:length(scobj.harmony.split)) {
  Singlet <- c(Singlet, 
               rownames(scobj.harmony.split[[i]]@meta.data) [scobj.harmony.split[[i]]@meta.data$DF.classifications_0.25 == "Singlet"])
}
finalcol <- ncol(scobj.harmony.split[[1]]@meta.data)
for (i in 1:length(scobj.harmony.split)) {
  all(scobj.harmony.split[[i]]@meta.data[rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet,   finalcol] == "Singlet") %>% stopifnot
  all(scobj.harmony.split[[i]]@meta.data[!(rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet),finalcol] != "Singlet") %>% stopifnot
}

scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)

scobj.harmony$dblfinder <- NA
scobj.harmony$dblfinder <- "doublet"
scobj.harmony$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"

dim(scobj.harmony)
# [1] 19041 16709
table(scobj.harmony$dblfinder)
# doublet singlet 
# 717   15992 
# Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1_scobj4run2(GSE193337).Rdata")

# 1000----
scobj.harmony.1000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  JoinLayers(assay = "RNA")

scobj.harmony.1000.20 <- scobj.harmony.1000 %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.1000.25 <- scobj.harmony.1000 %>%
  FindNeighbors(reduction = "pca", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.1000.30 <- scobj.harmony.1000 %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

pdf("run2_nf1000_h_noreg clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
dev.off()

save(scobj.harmony.1000, 
     scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
     file = "run2_nf1000_h_noreg scobj.harmony.Rdata")

# 1500----
scobj.harmony.1500 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  JoinLayers(assay = "RNA")

scobj.harmony.1500.20 <- scobj.harmony.1500 %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.1500.25 <- scobj.harmony.1500 %>%
  FindNeighbors(reduction = "pca", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.1500.30 <- scobj.harmony.1500 %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

pdf("run2_nf1500_h_noreg clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
dev.off()

save(scobj.harmony.1500, 
     scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
     file = "run2_nf1500_h_noreg scobj.harmony.Rdata")

# 2000----
scobj.harmony.2000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  JoinLayers(assay = "RNA")

scobj.harmony.2000.20 <- scobj.harmony.2000 %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.2000.25 <- scobj.harmony.2000 %>%
  FindNeighbors(reduction = "pca", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.2000.30 <- scobj.harmony.2000 %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

pdf("run2_nf2000_h_noreg clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
dev.off()

save(scobj.harmony.2000, 
     scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
     file = "run2_nf2000_h_noreg scobj.harmony.Rdata")

# select optimal parameter ----
scobj.harmony <- scobj.harmony.2000.30 %>% 
  RunUMAP(reduction = "pca", dims = 1:30)
# Annotation====================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.1")

epithelial <- c("EPCAM", "SFN", "KRT5", "KRT8", "KRT14","KLK3") # "SPRR3" 
endothelial <- c("VWF", "PECAM1", "ENG", "CDH5") #, "CCL14"
fibroblast <- c("FN1", "DCN", "COL1A1", "COL1A2") #, "COL3A1", "COL6A1"
SMC <- c("TAGLN", "CNN1", "PRKG1", "FOXP2")
pericyte <- c("RGS5", "MCAM", "ACTA2", "MYH11")
T_cell <- c("CD2", "CD3D", "CD3E", "CD3G") # 0 4 10 14 
B_cell <-	c("CD19", "CD79A", "MS4A1", "CD79B") # 12
plasma <-	c("JCHAIN", "MZB1", "IGHG1", "SDC1") # 8 16 , "CD79A"
myeloid <- c("CD68", "CD163", "LYZ", "CD14","FCGR3A","C1QA","C1QB","CST3") # , "IL3RA", "LAMP3", "CLEC4C", "GCA"
Neutrophil <- c("CSF3R","CXCL8","G0S2","IFITM2") # "FCGR3B","FPR1","BASP1","CXCR1","CXCR2","S100A11"
DC <- c('CCR7',	'CLEC9A', 'CD1C',	'IRF7',	'LILRA4')
Mast <-	c('TPSAB1','CPA3','HPGDS','KIT')# 'VWA5A','SLC18A2','HDC','CAPG','RGS13','IL1RL1','FOSB','GATA2'
Neural <- c("PLP1","NRNX1","NRNX2","NRNX3")
Proliferation <- c("TOP2A","BIRC5","MKI67","PCNA")
NK <- c('NKG7', 'GNLY', 'KLRD1', 'PRF1', 'FCGR3A', 'TYROBP', 'FCER1G')

p <- DotPlot(scobj.harmony, features = c(
  epithelial, endothelial, fibroblast, SMC, pericyte,
  T_cell, B_cell,plasma, myeloid, Neutrophil, Mast, DC, NK
) %>% unique, 
group.by = "RNA_snn_res.1") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = "Dotplot marker (GSE193337).pdf", width=9, height=8, units="in")

C2 <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.1", ident.1 = 2)
C18<- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.1", ident.1 = 18)
C20<- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.1", ident.1 = 20)
C23<- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.1", ident.1 = 23)
C25<- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.1", ident.1 = 25)

write.csv(C2, "C2.csv")
write.csv(C20, "C20.csv")
write.csv(C23, "C23.csv")
write.csv(C25, "C25.csv")

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(0)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(1)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(2)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(3)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(4)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(5)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(6)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(7)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(8)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(9)] <- "Mast"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(10)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(11)]<- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(12)]<- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(13)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(14)]<- "Mesenchymal"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(15)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(16)]<- "Mesenchymal"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(17)]<- "Myeloid" 
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(18)]<- "T cell" # NKT-like
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(19)]<- "Mesenchymal"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(20)]<- "Doublets"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(21)]<- "B cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(22)]<- "Mesenchymal"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(23)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(24)]<- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(25)]<- "NK"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(26)]<- "Endothelial"

save(scobj.harmony, file = "scobj.harmony(singlet_nf2000_h_noreg_pc30_res1)anno.Rdata")
load("scobj.harmony(singlet_nf2000_h_noreg_pc30_res1)anno.Rdata")

table(scobj.harmony$RNA_snn_res.1)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 1762 1639 1391 1231 1207 1010  805  776  667  593  572  485  465  430  400  375  357  307  293 
# 19   20   21   22   23   24   25   26 
# 219  215  167  161  150  119  112   84 
table(scobj.harmony$cell_type)
# B cell    Doublets Endothelial  Epithelial        Mast Mesenchymal     Myeloid          NK 
# 167         215        2096        4361         593        1137        1257         112 
# T cell 
# 6054 
# CNV --------------------------------------------------------------------------
setwd("copykat/")

test_sample = scobj.harmony$orig.ident[scobj.harmony$group == "Tumor"] %>% unique
copykat_list = lapply(test_sample, function(gsm){
  
  tmp  <- subset(
    scobj.harmony, 
    subset = orig.ident %in% gsm &
             cell_type  %in% c("Epithelial", "T cell", "B cell", "Myeloid")
    )
  exp.rawdata <- as.matrix(tmp@assays$RNA$raw_count)
  
  meta <- tmp@meta.data
  ref  <- rownames(meta)[meta$cell_type != "Epithelial"]
  
  copykat.test <- copykat(
    rawmat    = exp.rawdata, 
    id.type   = "S", 
    ngene.chr = 5, 
    win.size  = 25, 
    KS.cut    = 0.1, 
    sam.name  = gsm, 
    distance  = "euclidean", 
    norm.cell.names = ref,
    output.seg= FALSE, 
    plot.genes= FALSE, 
    genome    = "hg20",
    n.cores  = 1)
  
  meta$cell.names = rownames(meta)
  res = meta %>% 
    left_join(copykat.test$prediction, by = "cell.names")
  
  return(res)
})

save(copykat_list, file = "copykat_list.Rdata")

copykat_res  = copykat_list %>% bind_rows

meta.data <- scobj.harmony@meta.data
meta.data$cell.names <- rownames(meta.data) 
meta.data <- meta.data %>% 
  left_join(
    copykat_res %>% 
      dplyr::select(cell.names, copykat.pred), 
    by = "cell.names") %>% 
  as.data.frame
rownames(meta.data) <- meta.data$cell.names
scobj.harmony@meta.data <- meta.data

setwd("..")
save(scobj.harmony, file = "scobj.harmony_anno_copykat.Rdata")
# Figure2=======================================================================
fig = subset(scobj.harmony, subset = cell_type != "Doublets")
fig$copykat.pred[fig$cell_type != "Epithelial"] = NA
fig$cell_type2 = fig$cell_type
fig$cell_type2[fig$copykat.pred == "aneuploid"] = "Malignant"

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

pdf("run2_UMAP.pdf")
DimPlot(fig, reduction = "umap", group.by = "cell_type", label = TRUE, cols = colors_list)
DimPlot(fig, reduction = "umap", group.by = "cell_type2", label = TRUE, cols = colors_list)
DimPlot(fig, reduction = "umap", group.by = "orig.ident", label = TRUE)
DimPlot(fig, reduction = "umap", group.by = "sample", label = TRUE)
DimPlot(fig, reduction = "umap", group.by = "group", label = TRUE)
DimPlot(fig, reduction = "umap", group.by = "disease", label = TRUE)
DimPlot(fig, reduction = "umap", group.by = "copykat.pred", label = TRUE)
dev.off()

p_vln <- VlnPlot(fig, features = c("ACTB"), ncol = 4,
                 group.by = "cell_type2", cols = colors_list) + guides(col = "none")
ggsave(p_vln, filename = "VlnPlot.pdf", width=3, height=4, units="in")
# Pseudo bulk ==================================================================
Epi <- subset(fig, subset = cell_type == "Epithelial")
count_mat <- GetAssayData(Epi, assay = "RNA", layer = "counts")

meta <- Epi@meta.data

meta %>% 
  group_by(orig.ident) %>% 
  mutate(sum = n()) %>% 
  distinct(sum)

sample_info <- meta %>%
  select(orig.ident, group, disease, sample) %>%
  distinct(., .keep_all = TRUE)

pb_counts <- lapply(unique(Epi$orig.ident), function(orig.ident){
  barcodes = meta$cell.names[meta$orig.ident == orig.ident]
  pb_count = data.frame(
    V1 = rowSums(as.matrix(count_mat)[,barcodes]),
    row.names = rownames(count_mat)
  )
  colnames(pb_count) = orig.ident
  
  return(pb_count)
}) %>% bind_cols()


count <- pb_counts[,sample_info$orig.ident]

stopifnot(all(colnames(count) == sample_info$orig.ident))

condition = factor(sample_info$group, levels = c("Normal","Tumor"))
sample    = factor(sample_info$sample,levels = paste0("patient ", c(1:4)), labels = paste0("p", c(1:4)))
coldata   = data.frame(
  row.names = colnames(count), 
  condition = condition,
  sample    = sample)
dds       = DESeqDataSetFromMatrix(countData = count,
                              colData = coldata,
                              design = ~sample + condition)
dds$condition<- relevel(dds$condition, ref = "Normal") # 指定哪一组作为对照组
dds <- DESeq(dds)  
DEG <- results(dds, name="condition_Tumor_vs_Normal", independentFiltering = FALSE) %>%
  as.data.frame %>% 
  na.omit

save(DEG, file = "DEG.Rdata")

