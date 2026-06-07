rm(list = ls());gc()

setwd("")

# GSM5848795	Normal kidney allograft rep 1
# GSM5848796	Normal kidney allograft rep 2
# GSM5848797	Normal kidney allograft rep 3
# GSM5848798	Fibrotic kidney allograft rep 1
# GSM5848799	Fibrotic kidney allograft rep 2
# GSM5848800	Fibrotic kidney allograft rep 3
# GSM5848801	Fibrotic kidney allograft rep 4
# GSM5848802	Fibrotic kidney allograft rep 5
# GSM5848803	Fibrotic kidney allograft rep 6

# Load pkgs---------------------------------------------------------------------
set.seed(1011)
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
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

myqc4seurat <- function(seurat.obj,
                        xintercept1 = c(200,300,400,500,1000,5000,6000,7000,8000),
                        xintercept2 = c(200,500,1000,5000,10000,15000)){
  p1 <- VlnPlot(seurat.obj, features = c("nFeature_RNA"), group.by = "orig.ident", ncol = 1) + scale_y_continuous(breaks = c(200,500, 1000,2000,4000,6000,8000,10000))
  p2 <- VlnPlot(seurat.obj, features = c("nCount_RNA"), group.by = "orig.ident", ncol = 1)
  p3 <- VlnPlot(seurat.obj, features = c("percent.mt"), group.by = "orig.ident", ncol = 1) + scale_y_continuous(breaks = c(10,20))
  p4 <- VlnPlot(seurat.obj, features = c("percent.rp"), group.by = "orig.ident",  ncol = 1) + scale_y_continuous(breaks = c(10,20))
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
# Load data---------------------------------------------------------------------
set.seed(1011)
setwd("GSE195718_RAW/")
# options(future.globals.maxSize = 1000 * 1024^2) # 设置为1GB

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

hub <- dir()
names(hub) <- hub

gse195718 <- hub %>%
  Read10X %>%
  CreateSeuratObject(
    min.cells = 0,
    min.features = 0,
    project = "gse195718",
    assay = "RNA")

normal   <- c("GSM5848795", "GSM5848796", "GSM5848797")
fibrotic <- c("GSM5848798", "GSM5848799", "GSM5848800",
              "GSM5848801", "GSM5848802", "GSM5848803")
gse195718$group <- NA
gse195718$group[gse195718$orig.ident %in% normal]   <- "Normal"
gse195718$group[gse195718$orig.ident %in% fibrotic] <- "Fibrotic"
# I - run1 QC-------------------------------------------------------------------
scobj <- gse195718
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

pdf("run1 QC前.pdf")
myqc4seurat(seurat.obj = scobj)
dev.off()

scobj.qc <- subset(scobj, subset = orig.ident != "GSM5848799")
scobj.qc <- subset(
  scobj.qc, 
  subset = 
    nFeature_RNA > 500 &
    nFeature_RNA < 4000 & 
    nCount_RNA > 500 &
    nCount_RNA < 10000 &
    percent.mt < 10 &
    percent.rp < 10 &
    percent.hb < 1)  

pdf("run1 QC后.pdf")
myqc4seurat(seurat.obj = scobj.qc)
dev.off()

table(scobj.qc$orig.ident)
# GSM5848795 GSM5848796 GSM5848797 GSM5848798 GSM5848800 GSM5848801 GSM5848802 
# 7157       6589       3885       5971       5253       5242       3034 
# GSM5848803 
# 4245

dim(scobj.qc)
# [1] 36601 41376
# I - run1 raw data stat-----------------------------------------------------------------
gtf_data <- rtracklayer::import("refdata-gex-GRCh38-2020-A genes.gtf") %>% as.data.frame

keep <- gtf_data$gene_name[gtf_data$gene_type != 'lncRNA'] %>% unique
scobj<- scobj.qc[rownames(scobj.qc) %in% keep,]

# dim(scobj@assays$RNA$counts)
# [1] 20034 41376

table(scobj$orig.ident)
# GSM5848795 GSM5848796 GSM5848797 GSM5848798 GSM5848800 GSM5848801 GSM5848802 GSM5848803 
# 7157       6589       3885       5971       5253       5242       3034       4245 
# I - run1 blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(gse195718))
# [1] 239
scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj)
# [1] 19806 41376

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "run1 scobj.Rdata")
# I - run1 Nomalization & harmony & cluster----------------------------------------------
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

pdf("run1 nf2000_h_clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony, reduction = "pca", ndims = 50)
clustree(scobj.harmony.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.30, prefix = "RNA_snn_res.")
dev.off()

scobj.harmony <- scobj.harmony.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# I - run1 Visualization-----------------------------------------------------------------
pdf("run1 nf2000_h_pc30 umap.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.4", reduction = "umap", label = T)
FeaturePlot(scobj.harmony, features = "percent.mt")
dev.off()
# I - run1 Cell Cycle--------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.harmony <- CellCycleScoring(scobj.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scobj.harmony <- RunPCA(scobj.harmony, features = c(s.genes, g2m.genes), reduction.name = "cellcycle")

pdf("run1 nf2000_h_pc30 umap cellcycle.pdf")
RidgePlot(scobj.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
DimPlot(scobj.harmony, reduction = "cellcycle")
DimPlot(scobj.harmony, reduction = "umap")
dev.off()
# I - run1 doublets----------------------------------------------------------------------
save(scobj.harmony, file = "run1 nf2000_h_pc30 scobj.harmony.Rdata")

scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list

nPC = 30
for (i in 1:length(scobj.harmony.split)) {
  # pK Identification (ground-truth) -------------------------------------------
  sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:nPC)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
  ## Homotypic Doublet Proportion Estimate -------------------------------------
  homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]$RNA_snn_res.0.4)  
  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi <- round(0.1 * nrow(scobj.harmony.split[[i]]@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies ----------------
  # scobj.harmony.split[[i]] <- doubletFinder(scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/228
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/225#issuecomment-2786505997
  scobj.harmony.split[[i]] <- doubletFinder(
    scobj.harmony.split[[i]], PCs = 1:nPC, pN = 0.25, pK = pK, 
    nExp = nExp_poi.adj, sct = FALSE)
}

save(scobj.harmony.split, file = "run1 nf2000_h_pc30 scobj.harmony.split.Rdata")


Singlet <- c()
for (i in 1:length(scobj.harmony.split)) {
  Singlet <- c(Singlet, 
               rownames(scobj.harmony.split[[i]]@meta.data) [scobj.harmony.split[[i]]@meta.data$DF.classifications_0.25 == "Singlet"])
}
finalcol <- ncol(scobj.harmony.split[[1]]@meta.data)
for (i in 1:length(scobj.harmony.split)) {
  all(scobj.harmony.split[[i]]@meta.data[rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet,   finalcol] == "Singlet") %>% print
  all(scobj.harmony.split[[i]]@meta.data[!(rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet),finalcol] != "Singlet") %>% print
}

scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)

# scobj.h.dblfinder <- subset(scobj.harmony, subset = id %in% Singlet)
scobj.harmony$dblfinder <- NA
scobj.harmony$dblfinder <- "doublet"
scobj.harmony$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"

dim(scobj.harmony)
# [1] 19806 41376
table(scobj.harmony$dblfinder)
# doublet singlet 
# 3741   37635
# I - run2 Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1 scobj4run2.Rdata")

# I - run2 noreg========================================================================
# I - run2 1000 without regress cell cycle----
set.seed(1011)

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
save(scobj.harmony.1000,
     scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
     file = "run2 nf1000_h_noreg scobj.Rdata")

# I - run2 1500 without regress cell cycle----
set.seed(1011)

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
save(scobj.harmony.1500,
     scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
     file = "run2 nf1500_h_noreg scobj.Rdata")

# I - run2 2000 without regress cell cycle----
set.seed(1011)

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
save(scobj.harmony.2000,
     scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
     file = "run2 nf2000_h_noreg scobj.Rdata")


pdf("run2 nf1000_h_noreg_clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
dev.off()

pdf("run2 nf1500_h_noreg_clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
dev.off()

pdf("run2 nf2000_h_noreg_clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
dev.off()

# I - run2 Run Umap----
scobj.harmony <- scobj.harmony.1000.25 %>% 
  RunUMAP(reduction = "harmony", dims = 1:25)
# I - run2 Visualization=================================================================
pdf("run2 nf1000_h_pc25_noreg_res0.2 umap.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T) + NoLegend() 
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.2", reduction = "umap", label = T)
dev.off()
# I - run2 Annotation====================================================================
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
cell_cycle <- c("TOP2A","PCNA","MKI67","BIRC5")

if(TRUE){
  ####epithelial cells----
  renal_corpuscle <- c('PTPRQ', 'WT1', 'NTNG1', 'NPHS1', 'NPHS2', 'CLIC5', 'PODXL')
  glomerulus <- c('CLDN1', 'VCAM1', 'CFH', 'RBFOX1', 'ALDH1A2') 
  proximal_tubules <- c(
    'LRP2', 'CUBN', 'SLC13A1', 'SLC5A12', 'SLC13A3', 'SLC22A6', 'PRODH2', 'SLC5A2', 
    'SLC22A8', 'SLC34A1', 'SLC22A7', 'MOGAT1', 'SLC5A11', 'SLC22A24', 'SLC7A13', 
    'SLC5A8', 'ABCC3', 'SATB2') 
  PT <- c("MIOX", "ALDOB", "FABP1", "PCK1", "ANPEP")
  intermediate_tubules <- c(
    'CRYAB', 'TACSTD2', 'SLC44A5', 'KLRG2', 'COL26A1', 'BOC', 'VCAM1', 'SLC39A8', 
    'AQP1', 'LRRC4C', 'LRP2', 'UNC5D', 'SATB2', 'JAG1', 'ADGRL3', 'ID1', 'CLDN1', 
    'AKR1B1', 'CLDN4', 'BCL6', 'SH3GL3', 'SLC14A2', 'SMOC2', 'BCAS1', 'CLCNKA', 
    'CLDN10', 'PROX1') 
  Distal_tubules <- c(
    "CASR","SLC12A1","UMOD","NELL1","ESRRB","EGF","CLDN14","PROX1","MFSD4A","KCTD16",  
    "RAP1GAP","ANK2","CYFIP2","PPM1E","GP2","ENOX1","TMEM207","TMEM52B","CLDN16","WNK1M",   
    "NOS1","ROBO2","CALCR","PPFIA2","PAPPA2","SLC12A3","CNNM2","FGF13","KLHL3","LHX1",  
    "TRPM6","TRPM7","ADAMTS17","ITPKB","ZNF385D","HS6ST2","TRPV5","SLC8A1","SCN2A","HSD11B2", 
    "CALB1") %>% sort 
  Collecting_tubules <- c(
    'SLC8A1', 'SCN2A', 'HSD11B2', 'CALB1', 'KITLG', 'PCDH7', 'RALYL', 'TOX', 'SGPP1', 
    'SCNN1G', 'SCNN1B', 'KCNIP1', 'GATA3', 'AQP2', 'AQP3', 'FXYD4', 'SOX5', 'PDE10A', 
    'SLC25A29', 'ST6GAL1', 'PAPPA', 'SYK', 'FAM81A', 'PROM1', 'KCNK13', 'FXYD4', 'SOX5',
    'PHACTR1', 'PCDH7', 'SLC14A2', 'HS3ST5', 'TACSTD2', 'TP63', 'GPX2', 'FXYD3', 'KRT5',
    'ATP6V0D2', 'ATP6V1C2', 'TMEM213', 'CLNK', 'SLC4A1', 'SLC26A7', 'HS6ST3', 'NXPH2', 
    'LEF1', 'ADGRF5', 'CALB1', 'KIT', 'AQP6', 'STAP1', 'FAM184B', 
    'CALCA', 'SLC4A9', 'SLC35F3', 'SLC26A4', 'INSRR', 'TLDC2'
  ) %>% unique
  #### Endothelial Cell----
  #### Vascular Smooth Muscle Cell / Pericyte(stroma cells)----
  #### Fibroblast(stroma cells)----
  #### Immune Cells----
  # B Cell <- BANK1, BLK, MS4A1, BACH2
  # Plasma Cell <- IGKC, TENT5C, MZB1, FCRL5, CD38, JCHAIN
  # T Cell <- CD96, CD247, THEMIS, BCL11B, CAMK4, IL7R
  # Natural Killer T Cell <- CD96, CD247, RUNX3, GNLY, NKG7, CCL5, KLRF1, CCL4, GZMA
  # Mast Cell <- MS4A2, CPA3, KIT
  # M2 Macrophage <- F13A1, MRC1, CD163, STAB1, SLC1A3, CD14, FOLR2
  # Monocyte-derived Cell <- MSR1, ITGAX, HLA-DQA1, HLA_DRB1, CSF2RA, CD14, TRPM2
  # Classical Dendritic Cell <- ITGAX, HLA-DQA1, HLA-DRA, CSF2RA, CIITA, WDFY4, FLT3, ZNF366, CADM1, ZBTB46, CLEC9A
  # Plasmacytoid Dendritic Cell <- IRF8, CUX2, P2RY14, IL3RA, CLEC4C
  # Non-Classical Monocyte <- CTSS, IRAK3, TCF7L2, TNFRSF1B, FCN1, HLA-DRA, FCGR3A
  # Neutrophil <- S100A9, S100A8, IFITM2, FCGR3B, CD1C
}
markers <- readxl::read_excel("Renal_markers.xlsx")

marker_long <- markers %>%
  dplyr::select(
    subclass = `Subclass Level 1`,
    positive_markers = `Positive Markers`
  ) %>%
  filter(!is.na(positive_markers), positive_markers != "") %>%
  separate_rows(positive_markers, sep = ",") %>%
  mutate(
    gene = str_trim(positive_markers)
  ) %>%
  dplyr::select(subclass, gene) %>%
  filter(gene != "") %>%
  distinct()

head(marker_long)

p_dotplot <- DotPlot(scobj.harmony, features = c(
  epithelial, endothelial, fibroblast, SMC, pericyte,
  T_cell, B_cell,plasma, myeloid, Neutrophil, Mast, DC, cell_cycle,
  marker_long$gene
) %>% unique, 
group.by = "RNA_snn_res.0.5") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p_dotplot, filename = "DotPlot4Anno (res0.5) R2.pdf", width=8, height=40, units="in")

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(0)] <- "TAL" # Thick Ascending Limb of the loop of Henle
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(1)] <- "PT" # proximal tubules
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(2)] <- "DCT/CNT/PC"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(3)] <- "PT" # proximal tubules
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(4)] <- "TAL"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(5)] <- "PT"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(6)] <- "EC"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(7)] <- "TAL/DCT" # mixed
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(8)] <- "TAL"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(9)] <- "IC" # Collecting tubules
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(10)]<- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(11)]<- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(12)]<- "EC"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(13)]<- "FIB"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(14)]<- "IC" # Collecting tubules
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(15)]<- "POD" # renal corpuscle
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(16)]<- "VSM/P"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(17)]<- "PEC" # glomerulus
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(18)]<- "Plasma"

Clist <- list()
n <- 1
for (i in 0:18) {
  Clist[[n]] <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.5", ident.1 = as.character(i))
  write.csv(Clist[[n]], paste0("run2|nf1000_h_regccmt_pc25_res0.5_C",n-1,".csv"))
  
  n <- n + 1
}

# save(scobj.harmony, file = "scobj.harmony(singlet_nf1000_h_noreg_pc25_res0.2)anno.Rdata")
load("scobj.harmony(singlet_nf1000_h_noreg_pc25_res0.2)anno.Rdata")

table(scobj.harmony$RNA_snn_res.0.5)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
# 5501 4491 4172 3458 2621 2569 2568 2389 1761 1629 1335 1126  948  819  677  574  433  378 
# 18 
# 186

table(scobj.harmony@meta.data$cell_type)
# DCT/CNT/PC         EC        FIB         IC    Myeloid        PEC     Plasma        POD 
# 4172       3516        819       2306       1126        378        186        574 
# PT     T cell        TAL    TAL/DCT      VSM/P 
# 10518       1335       9883       2389        433 
# I - Visualization=============================================================
fig <- scobj.harmony
fig$sample <- fig$orig.ident
fig$RNA_snn_res.0.5 <- factor(fig$RNA_snn_res.0.5, levels = c(0:18))
# I - Visualization|umap----
p1 <- DimPlot(fig, reduction = "umap", group.by = "cell_type", cols = colors_list,
              label = TRUE, pt.size = 0.1, raster = FALSE)
p2 <- DimPlot(fig, reduction = "umap", group.by = "orig.ident",
              label = TRUE, pt.size = 0.1, raster = FALSE)
p3 <- DimPlot(fig, reduction = "umap", group.by = "group", cols = c('#4DBBD5', '#E64B35'),
              label = TRUE, pt.size = 0.1, raster = FALSE)
ggsave(p1, filename = "Figure1 UMAP(cell_type).pdf",  width=8.5, height=8, units="in")
ggsave(p2, filename = "Figure1 UMAP(orig.ident).pdf", width=8.5, height=8, units="in")
ggsave(p3, filename = "Figure1 UMAP(group).pdf",   width=8.5, height=8, units="in")
# I - Visualization|marker ----
p <- DotPlot(fig, features = c(
  endothelial, Mast, fibroblast, pericyte, SMC, myeloid, Neutrophil, T_cell,
  marker_long$gene) %>% unique, 
  group.by = "cell_type") + 
  scale_color_viridis() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = "Figure1 MARKER.pdf", width=8, height=40, units="in")
# I - Visualization|sankey by sample----
dat <- fig@meta.data %>% 
  dplyr::select(group, cell_type) %>% 
  mutate(group = factor(group, levels = c("Normal","Fibrotic")))
p <- dat %>% 
  mutate(value = 1) %>% 
  group_by(group, cell_type) %>%
  summarise(n = n()) %>% 
  group_by(group) %>%
  mutate(sum = sum(n)) %>% 
  mutate(n = n / sum) %>% 
  ggplot(aes(x = group, y = n, fill = cell_type, stratum = cell_type, alluvium = cell_type))+
  geom_stratum(width = 0.5, color='white') +
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear") +
  scale_color_manual(values = colors_list) +
  scale_fill_manual(values = colors_list) +
  labs(x = "", y = "Percent", fill = "Cell type") +
  # coord_flip() +
  theme_bw()+
  theme(text = element_text(size = 20))
ggsave(p, filename = "Figure1 percent.pdf",
       width=6, height=6, units="in")
# I - Visualization|sankey by group ----
dat <- fig@meta.data %>% 
  dplyr::select(sample, group, cell_type) %>% 
  mutate(group = factor(group, levels = c("Normal","Fibrotic")))
p <- dat %>% 
  mutate(value = 1) %>% 
  group_by(sample, group, cell_type) %>%
  summarise(n = n()) %>% 
  group_by(sample, group) %>%
  mutate(sum = sum(n)) %>% 
  mutate(n = n / sum) %>% 
  ggplot(aes(x = sample, y = n, fill = cell_type, stratum = cell_type, alluvium = cell_type))+
  geom_stratum(width = 0.5, color='white') +
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear") +
  scale_color_manual(values = colors_list) +
  scale_fill_manual(values = colors_list) +
  labs(x = "", y = "Percent", fill = "Cell type") +
  # coord_flip() +
  theme_bw()+
  theme(text = element_text(size = 20))
ggsave(p, filename = "Figure1 percent sample.pdf",
       width=8, height=6, units="in")
# I - Visualization|stat ----
table(fig$sample, fig$cell_type) %>%
  as.data.frame %>%
  tidyr::pivot_wider(names_from = "Var2", values_from = "Freq") %>%
  write.csv("Figure1(number).csv")

table(fig$sample, fig$cell_type) %>%
  as.data.frame %>%
  group_by(Var1) %>%
  mutate(sum = sum(Freq)) %>% 
  mutate(prop= Freq / sum) %>% 
  dplyr::select(-Freq, -sum) %>% 
  tidyr::pivot_wider(names_from = "Var2", values_from = "prop") %>%
  write.csv("Figure1(prop).csv")

data <- fig@meta.data
data$majorCluster = data$cell_type
data$patient = data$sample
data$loc = data$group
Roe <- calTissueDist(data,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
# Roe
# #                     APE       CTEPH
# # Endothelial 0.048567665 1.851972463
# # Mast        0.007868925 1.888416679
# # Mesenchyme  0.086305570 1.818179565
# # Myeloid     1.207129289 0.814523384
# # Neutrophil  1.607151724 0.456318090
# # T cell      0.089134533 1.815646333

df_long <- Roe %>%
  as.data.frame %>%
  transmute(Cell_Type = Var1,
            Disease_Type = Var2,
            Roe = Freq)
p <- ggplot(df_long, aes(x = Disease_Type, y = Cell_Type, fill = Roe)) +
  geom_tile(color = "white", size = 0.5) +  # 白色格子边框让热图更清晰
  geom_text(aes(label = sprintf("%.3f", Roe)), color = "black", size = 4) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", 
                       midpoint = 1.0, # 观察数据，中位数设在 1.0 左右最适合展现两极分化
                       name = "Value") +
  labs(x = "", y = "", title = "Roe") +
  theme_minimal() + # 使用干净的极简主题
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 11, color = "black", face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank() # 移除网格线
  )
ggsave(p, filename = "Figure1 Roe.pdf", width=5, height=5, units="in")
# I - DEA|pseudobulk------------------------------------------------------------
normal   <- c("GSM5848795", "GSM5848796", "GSM5848797")
fibrotic <- c("GSM5848798", "GSM5848799", "GSM5848800",
              "GSM5848801", "GSM5848802", "GSM5848803")
tb <- table(fig$cell_type)
setwd("DEA/")
for (i in names(tb)[tb >= 100]) {
  tmp      <- fig %>% subset(subset = cell_type == i)
  
  counts   <- AggregateExpression(tmp, assays = "RNA", group.by = "orig.ident")
  counts   <- counts$RNA %>% as.data.frame
  
  group    <- ifelse(colnames(counts) %in% normal, "Normal", 
                     ifelse(colnames(counts) %in% fibrotic,"Fibrotic", NA))
  if(sum(group == "Normal") == 0 | sum(group == "Fibrotic") == 0) {
    message(paste0("Not enough sample in ", i, " cell!"))
    next
  }
  
  keep     <- rowSums(counts >= 10) >= 2
  counts   <- counts[keep, ]
  

  condition<- factor(group, levels = c("Normal","Fibrotic"))
  coldata  <- data.frame(row.names = colnames(counts), condition)
  dds      <- DESeqDataSetFromMatrix(
    countData = counts %>% as.matrix,
    colData = coldata,
    design = ~condition)
  dds      <- DESeq(dds)  
  DEG      <- results(dds, cooksCutoff=FALSE,
                 name="condition_Fibrotic_vs_Normal",
                 independentFiltering = TRUE) %>% 
    as.data.frame %>% 
    arrange(padj) %>% 
    na.omit
  
  i <- str_replace_all(i, "/", "_")
  
  save(counts, file = paste0("pseudobulk_",i,".Rdata"))
  save(DEG, file = paste0(i,"_DESeq2-DEGs.Rdata"))
  export::table2excel(DEG, add.rownames = TRUE, sheetName = i, paste0(i,"_DESeq2-DEGs.xlsx"))
  
  print(i)
}

pList <- list()
n <- 1
for (i in c("B cell","Endothelial","Epithelial","Mast","Mesenchyme","Myeloid",
            "Neutrophil","Plasma","T cell")) {
  load(paste0("1/1.DEA/bulk 拆分SMC/drinking/",i,"_DESeq2-DEGs.Rdata"))
  p <- geom_volcano(DEG_drinker_vs_nondrinker, title = i)
  pList[[n]] <- p
  n <- n + 1
}

p_out <- patchwork::wrap_plots(pList) +
  plot_layout(ncol = 3, nrow = 3)
ggsave(p_out, filename = paste0(i,"volcanoplot drinking.pdf"), width=15, height=15, units="in")
# I - Visualization|Featureplot-------------------------------------------------------------------
p_out <- lapply(paste0("GRK", c(1:7)), 
                function(x){FeaturePlot(fig, features = x, pt.size = 0.01) + 
                    scale_colour_gradientn(
                      colours = colorRampPalette(
                        c('#eeeeee11','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410',
                          '#E52C22','#DD1D23','#C20030','#930039','#8C003A',
                          '#6F003D','#56033F'))(1000))}) %>% 
  patchwork::wrap_plots(ncol = 4, nrow = 2)
ggsave(p_out, filename = "FeaturePlot GRK6.pdf", width=12, height=5, units="in")

p_out <- VlnPlot(fig, features = paste0("GRK", c(1:7)), group.by = "cell_type", cols = colors_list)
ggsave(p_out, filename = "VlnPlot GRK6.pdf", width=10, height=10, units="in")
