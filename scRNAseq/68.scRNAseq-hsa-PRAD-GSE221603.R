GSM6890190	mHSPC_02
GSM6890192	mHSPC_03
GSM6890194	Pca_01_N
GSM6890196	Pca_01_T1
GSM6890198	Pca_01_T2
GSM6890200	Pca_02_T
GSM6890202	Pca_03_TA
GSM6890204	Pca_03_TB
GSM6890206	Pca_04_N
GSM6890208	Pca_04_T1
GSM6890209	Pca_06_N
GSM6890210	Pca_06_T

mHSPC_02	GSM6890190	mHSPC	   batch3	 mHSPC_02
mHSPC_03	GSM6890192	mHSPC	   batch3	 mHSPC_03
Pca_01_N	GSM6890194	Normal   batch1	 Pca_01
Pca_01_T1	GSM6890196	Tumor	   batch2	 Pca_01
Pca_01_T2	GSM6890198	Tumor	   batch2	 Pca_01
Pca_02_T	GSM6890200	Tumor	   batch1	 Pca_02
Pca_03_TA	GSM6890202	Tumor	   batch2	 Pca_03
Pca_03_TB	GSM6890204	Tumor	   batch2	 Pca_03
Pca_04_N	GSM6890206	Normal   batch1	 Pca_04
Pca_04_T1	GSM6890208	Tumor	   batch2	 Pca_04
Pca_06_N	GSM6890209	Normal   batch1	 Pca_06
Pca_06_T	GSM6890210	Tumor	   batch2	 Pca_06


setwd("GSE221603/")

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
colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

hub <- dir()
names(hub) <- hub

scobj <- hub %>%
  Read10X %>%
  CreateSeuratObject(
    min.cells = 0,
    min.features = 0,
    project = "GSE221603",
    assay = "RNA")

scobj
# An object of class Seurat 
# 36601 features across 34112 samples within 1 assay 
# Active assay: RNA (36601 features, 0 variable features)
# 1 layer present: counts

scobj <- scobj[!str_detect(rownames(scobj), "^A[CFJLP][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^B[X][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^ENSG"),]

metadata <- scobj@meta.data
metadata <- metadata %>% 
  mutate(
    sample = case_when(orig.ident %in% c("GSM6890194","GSM6890196","GSM6890198") ~ "P1",
                       orig.ident %in% c("GSM6890200") ~ "P2",
                       orig.ident %in% c("GSM6890202","GSM6890204") ~ "P3",
                       orig.ident %in% c("GSM6890206","GSM6890208") ~ "P4",
                       orig.ident %in% c("GSM6890209","GSM6890210") ~ "P6",
                       orig.ident %in% c("GSM6890190") ~ "mHSPC_02",
                       orig.ident %in% c("GSM6890192") ~ "mHSPC_03"),
    tissue = case_when(orig.ident %in% c("GSM6890194","GSM6890206","GSM6890209") ~ "NAT", # Normal_Adjacent_Tissue
                       orig.ident %in% c("GSM6890190", "GSM6890192") ~ "mHSPC", # metastatic hormone sensitive prostate cancer
                       TRUE ~ "Tumor"),
    batch  = case_when(orig.ident %in% c("GSM6890194","GSM6890200","GSM6890206","GSM6890209") ~ "batch1",
                       orig.ident %in% c("GSM6890196","GSM6890198","GSM6890202","GSM6890204","GSM6890208","GSM6890210") ~ "batch2",
                       orig.ident %in% c("GSM6890190", "GSM6890192") ~ "batch3"))
metadata -> scobj@meta.data

dim(scobj@assays$RNA$counts)
# [1] 24249 34112

table(scobj@meta.data$orig.ident)
# GSM6890190 GSM6890192 GSM6890194 GSM6890196 GSM6890198 GSM6890200 GSM6890202 GSM6890204 
# 6011       1884       1566       1895        567       3570       4646       5535 
# GSM6890206 GSM6890208 GSM6890209 GSM6890210 
# 1337       1998       1915       3188 
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
    nFeature_RNA < 6000 &
    nCount_RNA > 500 &
    nCount_RNA < 10000 &
    percent.mt < 20 &
    percent.rp < 25 &
    percent.hb < 0.2)  

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
# GSM6890190 GSM6890192 GSM6890194 GSM6890196 GSM6890198 GSM6890200 GSM6890202 GSM6890204 
# 2945        879        555        507        266       2308        946       1089 
# GSM6890206 GSM6890208 GSM6890209 GSM6890210 
# 721       1460       1426        812 

dim(scobj)
# [1] 24249 13914
# blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(scobj))
# [1] 239

scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj)
# [1] 24011 13914

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "run1=scobj(GSE221603).Rdata")
# Nomalization & harmony & cluster----------------------------------------------
scobj.harmony <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "sample", "batch"),
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

pdf("run1=nf2000_h clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony, reduction = "pca", ndims = 50)
clustree(scobj.harmony.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.30, prefix = "RNA_snn_res.")
dev.off()

scobj.harmony <- scobj.harmony.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Visualization-----------------------------------------------------------------
pdf("run1=nf2000_h_pc30.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "sample", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "tissue", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "batch", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.7", reduction = "umap", label = T)
FeaturePlot(scobj.harmony, features = "percent.mt")
dev.off()
# Cell Cycle--------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.harmony <- CellCycleScoring(scobj.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scobj.harmony <- RunPCA(scobj.harmony, features = c(s.genes, g2m.genes), reduction.name = "cellcycle")

pdf("run1=cell cycle.pdf")
RidgePlot(scobj.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
DimPlot(scobj.harmony, reduction = "cellcycle")
DimPlot(scobj.harmony, reduction = "umap")
dev.off()
# doublets----------------------------------------------------------------------
save(scobj.harmony, file = "run1=scobj.harmony(GSE221603).Rdata")


require(DoubletFinder)
scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list

for (i in 1:length(scobj.harmony.split)) {
  # pK Identification (ground-truth) -------------------------------------------
  sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:30)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
  ## Homotypic Doublet Proportion Estimate -------------------------------------
  homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]@meta.data$RNA_snn_res.0.7)  
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

save(scobj.harmony.split, file = "run1=scobj.harmony.split(GSE221603).Rdata")

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
  rownames(scobj.harmony.split[[12]]@meta.data)[scobj.harmony.split[[12]]@meta.data$DF.classifications_0.25 == "Singlet"])

scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)

# scobj.h.dblfinder <- subset(scobj.harmony, subset = id %in% Singlet)
scobj.harmony@meta.data$dblfinder <- NA
scobj.harmony@meta.data$dblfinder <- "doublet"
scobj.harmony@meta.data$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"

dim(scobj.harmony)
# [1] 24011 13914
table(scobj.harmony@meta.data$dblfinder)
# doublet singlet 
# 596   13318
# Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1=scobj4run2(GSE221603).Rdata")

# 1000----
scobj.harmony.1000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData(vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "sample", "batch"),
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
pdf("run2=nf1000_h_regccmt clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony.1000, 
     scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
     file = "run2=nf1000_h_regccmt scobj.harmony.Rdata")

# 1500----
scobj.harmony.1500 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData(vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "sample", "batch"),
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
pdf("run2=nf1500_h_regccmt clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony.1500, 
     scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
     file = "run2=nf1500_h_regccmt scobj.harmony.Rdata")

# 2000----
scobj.harmony.2000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "sample", "batch"),
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
pdf("run2=nf2000_h_regccmt clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony.2000, 
     scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
     file = "run2=nf2000_h_regccmt scobj.harmony.Rdata")


scobj.harmony <- scobj.harmony.2000.20 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20)
# Annotation====================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.1")

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


p <- DotPlot(scobj.harmony, features = c(
  epithelial, endothelial, fibroblast, SMC, pericyte,
  T_cell, B_cell,plasma, myeloid, Neutrophil, Mast, DC, Neural,
  Proliferation
) %>% unique, 
group.by = "RNA_snn_res.1") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = "Dotplot marker (GSE221603).pdf", width=9, height=8, units="in")

C17 <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.1", ident.1 = 17)

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(0)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(1)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(2)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(3)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(4)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(5)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(6)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(7)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(8)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(9)] <- "Mysenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(10)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(11)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(12)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(13)]<- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(14)]<- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(15)]<- "Mast"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(16)]<- "Mysenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(17)]<- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(18)]<- "B cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(19)]<- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(20)]<- "Proliferation"

# save(scobj.harmony, file = "scobj.harmony(singlet_nf2000_h_regccmt_pc20_res1)anno.Rdata")
load("scobj.harmony(singlet_nf2000_h_regccmt_pc20_res1)anno.Rdata")

table(scobj.harmony$RNA_snn_res.1)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
# 2414 1870 1087  884  878  877  653  623  593  561  487  486  387  386  218  201  198  179 
# 18   19   20 
# 150   95   91 
table(scobj.harmony$cell_type)
# B cell   Endothelial    Epithelial          Mast       Myeloid    Mysenchyme 
# 150          1263          5119           201          1197           759 
# Proliferation        T cell 
# 91          4538 
# Figure2=======================================================================
require(patchwork)
p1 <- DimPlot(scobj.harmony, reduction = "umap", 
              group.by = "cell_type", 
              label = TRUE)
p2 <- DimPlot(scobj.harmony, reduction = "umap", 
              group.by = "RNA_snn_res.1", 
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




