set.seed(1011)
rm(list = ls());gc()

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')
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
# library(slingshot) # Trajectory
# library(tradeSeq) # Trajectory
# library(monocle) # Trajectory
library(monocle3) # Trajectory
library("CellChat") # Communication
# library("iTALK") # Communication
# library(infercnv) # CNV
# library(copykat) # CNV
# library(ArchR)
# Load data---------------------------------------------------------------------
hub <- c("Sample_01_H", "Sample_01_T", "Sample_02_H", 
         "Sample_02_T", "Sample_03_H", "Sample_03_T")
names(hub) <- hub

scobj <- hub %>%
  Read10X %>%
  CreateSeuratObject(
    min.cells = 0,
    min.features = 0,
    project = "GSE153892",
    assay = "RNA")

dim(scobj@assays$RNA$counts)
# [1]  38606 146635

scobj <- scobj[!str_detect(rownames(scobj), "^A[CFJLP][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^B[X][0-9]+.[0-9]+$"),]

dim(scobj@assays$RNA$counts)
# [1] 23562 17997

scobj@meta.data$orig.ident <- str_remove_all(rownames(scobj@meta.data), "_[ACTG]+-[0-9]$")
scobj@meta.data$sample <- str_remove_all(rownames(scobj@meta.data), "_[HT]_[ACTG]+-[0-9]$")
scobj@meta.data$tissue <- rownames(scobj@meta.data) %>% 
  str_remove_all("^Sample_[0-3]{2}_") %>% 
  str_remove_all("_[ACTG]+-[0-9]$")

table(scobj@meta.data$orig.ident)
# Sample_01_H Sample_01_T Sample_02_H Sample_02_T Sample_03_H Sample_03_T 
# 2598        3818        1222        4165        2560        3634 
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

pdf("QC.pdf")
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA"), ncol = 1)
VlnPlot(scobj, group.by = "orig.ident", features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.rp"), ncol = 1)
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 200 &
    nFeature_RNA < 2000 & 
    nCount_RNA > 200 &
    nCount_RNA < 4000 &
    percent.mt < 10 &
    percent.rp < 40 &
    percent.hb < 0.2)  

pdf("QC2.pdf")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(scobj, features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, features = c("percent.rp"), ncol = 1)
VlnPlot(scobj, features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

table(scobj@meta.data$orig.ident)
# Sample_01_H Sample_01_T Sample_02_H Sample_02_T Sample_03_H Sample_03_T 
# 2265        2909         969        3961        1805        2381 

# dim(scobj)
# [1]  25644 112871
# blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("C:/D/R project/Multi-Omics-Routine/scRNAseq/blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(scobj))
# [1] 243
# scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
# dim(scobj)
# [1]  25414 112871

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "scobj(GSE153892).Rdata")
# Nomalization & harmony & cluster----------------------------------------------
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

scobj.harmony.20 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.25 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:25) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))
scobj.harmony.30 <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

pdf("vst.2000.h clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony, reduction = "pca", ndims = 50)
clustree(scobj.harmony.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.30, prefix = "RNA_snn_res.")
dev.off()

scobj.harmony <- scobj.harmony.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Visualization-----------------------------------------------------------------
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.4")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.4", reduction = "umap", label = T)
# Cell Cycle--------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.harmony <- CellCycleScoring(scobj.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(scobj.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

scobj.harmony <- RunPCA(scobj.harmony, features = c(s.genes, g2m.genes), reduction.name = "cellcycle")
DimPlot(scobj.harmony, reduction = "cellcycle")
DimPlot(scobj.harmony, reduction = "umap")
# doublets----------------------------------------------------------------------
save(scobj.harmony, file = "scobj.harmony(GSE153892).Rdata")

if(TRUE){
  require(DoubletFinder)
  
  scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list
  
  for (i in 1:length(scobj.harmony.split)) {
    # pK Identification (ground-truth) -------------------------------------------
    sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:30)
    sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
    ## Homotypic Doublet Proportion Estimate -------------------------------------
    homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]@meta.data$RNA_snn_res.0.4)  
    ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi <- round(0.05 * nrow(scobj.harmony.split[[i]]@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    ## Run DoubletFinder with varying classification stringencies ----------------
    # scobj.harmony.split[[i]] <- doubletFinder(scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    
    # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/228
    # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/225#issuecomment-2786505997
    scobj.harmony.split[[i]] <- doubletFinder(
      scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, 
      nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  }
  
  save(scobj.harmony.split, file = "scobj.harmony.split(GSE153892).Rdata")
  
  Singlet <- c(
    rownames(scobj.harmony.split[[1]]@meta.data) [scobj.harmony.split[[1]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[2]]@meta.data) [scobj.harmony.split[[2]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[3]]@meta.data) [scobj.harmony.split[[3]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[4]]@meta.data) [scobj.harmony.split[[4]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[5]]@meta.data) [scobj.harmony.split[[5]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[6]]@meta.data) [scobj.harmony.split[[6]]@meta.data$DF.classifications_0.25 == "Singlet"])
  
  scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)
  
  # scobj.h.dblfinder <- subset(scobj.harmony, subset = id %in% Singlet)
  scobj.harmony@meta.data$dblfinder <- NA
  scobj.harmony@meta.data$dblfinder <- "doublet"
  scobj.harmony@meta.data$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"
  
  dim(scobj.harmony)
  # [1] 23562 14290
  table(scobj.harmony@meta.data$dblfinder)
  # doublet singlet 
  # 567   13723
}
# Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="scobj_for_2nd(GSE153892).Rdata")

# 1000 without regress cell cycle
if(TRUE){
  scobj.harmony.1000 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
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
  
  pdf("vst.1000.h_clustree 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.1000,
       scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
       file = "scobj.harmony.1000(without regression).Rdata")
}
# 1500 without regress cell cycle
if(TRUE){
  scobj.harmony.1500 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
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
  
  pdf("vst.1500.h_clustree 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.1500,
       scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
       file = "scobj.harmony.1500(without regression).Rdata")
}
# 2000 without regress cell cycle
if(TRUE){
  scobj.harmony.2000 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
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
  
  pdf("vst.2000.h_clustree 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.2000,
       scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
       file = "scobj.harmony.2000(without regression).Rdata")
}

scobj.harmony <- scobj.harmony.1000.25 %>% 
  RunUMAP(reduction = "harmony", dims = 1:25)
# Visualization=================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.5")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.5", reduction = "umap", label = T)
# Annotation====================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.5")

epithelial <- c("EPCAM", "SFN", "KRT5", "KRT14") #, "SPRR3"
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

DotPlot(scobj.harmony, features = c(
  T_cell, B_cell,plasma, myeloid, Neutrophil, Mast , DC
) %>% unique, 
group.by = "RNA_snn_res.0.5") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(0)] <- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(1)] <- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(2)] <- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(3)] <- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(4)] <- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(5)] <- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(6)] <- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(7)] <- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(8)] <-  
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(9)] <- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(10)]<- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(11)]<- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(12)]<- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(13)]<- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(14)]<- 
  scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(15)]<- 
  
  save(scobj.harmony, file = "scobj.harmony[1000-25-res0.5] without annot.Rdata")

table(scobj.harmony$RNA_snn_res.0.5)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 2507 2023 1683 1487 1399 1062  784  655  560  528  342  248  190  182   37   36 
table(scobj.harmony@meta.data$cell_type)

