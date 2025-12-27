GSM5793824	scRNA-seq of patient 1 benign sample
GSM5793825	scRNA-seq of patient 2 benign sample
GSM5793826	scRNA-seq of patient 3 benign sample
GSM5793827	scRNA-seq of patient 4 benign sample
GSM5793828	scRNA-seq of patient 1 tumor sample
GSM5793829	scRNA-seq of patient 2 tumor sample
GSM5793831	scRNA-seq of patient 3 tumor sample
GSM5793832	scRNA-seq of patient 4 tumor sample

setwd("GSE193337/")

rm(list = ls());gc()
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
    project = "GSE193337",
    assay = "RNA")

scobj
# An object of class Seurat 
# 33538 features across 39516 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
# 1 layer present: counts

scobj <- scobj[!str_detect(rownames(scobj), "^A[CFJLP][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^B[X][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^ENSG"),]

metadata <- scobj@meta.data
metadata <- metadata %>% 
  mutate(
    tissue = str_extract(rownames(metadata), "P[1-4]{1}[nt]{1}") %>% 
      str_remove_all("^P[1-4]{1}"))
metadata -> scobj@meta.data

dim(scobj@assays$RNA$counts)
# [1] 23562 39516

table(scobj@meta.data$orig.ident)
# GSM5793824 GSM5793825 GSM5793826 GSM5793827 GSM5793828 GSM5793829 GSM5793831 GSM5793832 
# 6674       2644       1875       6557       2574       5922       5358       7912
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
    nFeature_RNA < 5000 &
    nCount_RNA > 500 &
    nCount_RNA < 15000 &
    percent.mt < 20 &
    percent.rp < 30 &
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
# GSM5793824 GSM5793825 GSM5793826 GSM5793827 GSM5793828 GSM5793829 GSM5793831 GSM5793832 
# 5290       1045        767       3183       1549       2388       1625       2969

dim(scobj)
# [1] 23562 18810
# blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(scobj))
# [1] 243

scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj)
# [1] 23320 18810

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "run1=scobj(GSE193337).Rdata")
# Nomalization & harmony & cluster----------------------------------------------
scobj.harmony <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident"),
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
pdf("run1=nf2000_h_pc30 umap.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "tissue", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.8", reduction = "umap", label = T)
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
save(scobj.harmony, file = "run1=nf2000_h_pc30 scobj.harmony(GSE193337).Rdata")

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
    homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]@meta.data$RNA_snn_res.0.8)  
    ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi <- round(0.075 * nrow(scobj.harmony.split[[i]]@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    ## Run DoubletFinder with varying classification stringencies ----------------
    # scobj.harmony.split[[i]] <- doubletFinder(scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    
    # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/228
    # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/225#issuecomment-2786505997
    scobj.harmony.split[[i]] <- doubletFinder(
      scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, 
      nExp = nExp_poi.adj, sct = FALSE)
  }
  
  save(scobj.harmony.split, file = "run1=nf2000_h_pc30 scobj.harmony.split(GSE193337).Rdata")
  
  Singlet <- c(
    rownames(scobj.harmony.split[[1]]@meta.data) [scobj.harmony.split[[1]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[2]]@meta.data) [scobj.harmony.split[[2]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[3]]@meta.data) [scobj.harmony.split[[3]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[4]]@meta.data) [scobj.harmony.split[[4]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[5]]@meta.data) [scobj.harmony.split[[5]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[6]]@meta.data) [scobj.harmony.split[[6]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[7]]@meta.data) [scobj.harmony.split[[7]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[8]]@meta.data) [scobj.harmony.split[[8]]@meta.data$DF.classifications_0.25 == "Singlet"])
  
  scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)
  
  # scobj.h.dblfinder <- subset(scobj.harmony, subset = id %in% Singlet)
  scobj.harmony@meta.data$dblfinder <- NA
  scobj.harmony@meta.data$dblfinder <- "doublet"
  scobj.harmony@meta.data$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"
  
  dim(scobj.harmony)
  # [1] 23320 18810
  table(scobj.harmony@meta.data$dblfinder)
  # doublet singlet 
  # 1253   17557
}
# Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1=scobj4run2(GSE193337).Rdata")

# 1000
scobj.harmony.1000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData(vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident"),
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

# 1500
scobj.harmony.1500 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData(vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident"),
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

# 2000
scobj.harmony.2000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident"),
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

scobj.harmony <- scobj.harmony.1000.20 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20)
# Visualization=================================================================
pdf("run2=nf1000_h_pc res0.6.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "tissue", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.6", reduction = "umap", label = T)
dev.off()
# Annotation====================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.6")

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
Proliferation <- c("TOP2A","MKI67","PCNA","BIRC5","DLGAP5")

DotPlot(scobj.harmony, features = c(
  epithelial, endothelial, fibroblast, SMC, pericyte,
  T_cell, B_cell,plasma, myeloid, Neutrophil, Mast, DC, Neural,
  Proliferation
) %>% unique, 
group.by = "RNA_snn_res.0.6") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(0)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(1)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(2)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(3)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(4)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(5)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(6)] <- "Mysenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(7)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(8)] <- "Mast"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(9)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(10)]<- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(11)]<- "B cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(12)]<- "Proliferation"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(13)]<- "Mysenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.6 %in% c(14)]<- "Endothelial"

# save(scobj.harmony, file = "scobj.harmony(singlet_nf1000_h_pc20_regccmt_res_0.6)anno.Rdata")
load("scobj.harmony(singlet_nf1000_h_pc20_regccmt_res_0.6)anno.Rdata")

table(scobj.harmony$RNA_snn_res.0.6)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 3903 2862 1759 1692 1662 1487  822  681  628  620  596  388  204  192   61 
table(scobj.harmony$cell_type)
# B cell   Endothelial    Epithelial          Mast       Myeloid    Mysenchyme 
# 388          1820          4841           628          1277          1014 
# Proliferation        T cell 
# 204          7385 
# Figure2=======================================================================
require(patchwork)
p1 <- DimPlot(scobj.harmony, reduction = "umap", 
              group.by = "cell_type", 
              label = TRUE)
p2 <- DimPlot(scobj.harmony, reduction = "umap", 
              group.by = "orig.ident", 
              label = TRUE)
p_out1 <- p1 + p2
ggsave(p_out1, filename = "UMAP(GSE193337).pdf", width=10, height=5, units="in")

p_fea <- lapply(c(gene_symbols), 
                function(x){FeaturePlot(scobj.harmony, features = x) + 
                    scale_colour_gradientn(
                      colours = colorRampPalette(
                        c('gray90','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410',
                          '#E52C22','#DD1D23','#C20030','#930039','#8C003A',
                          '#6F003D','#56033F'))(1000))}) %>% 
  patchwork::wrap_plots(ncol = 2, nrow = 2)
p_vln <- VlnPlot(scobj.harmony, c(gene_symbols), group.by = "cell_type")
p_tissue <- VlnPlot(scobj.harmony, c(gene_symbols), group.by = "tissue")

ggsave(p1, filename = "UMAP(GSE193337).pdf", width=6, height=5, units="in")
ggsave(p_fea, filename = "FeaturePlot(GSE193337).pdf", width=8, height=8, units="in")
ggsave(p_vln, filename = "Vln_plot(GSE193337).pdf", width=9, height=6, units="in")
ggsave(p_tissue, filename = "Vln_tissue(GSE193337).pdf", width=9, height=6, units="in")

# sc.epithelial <- subset(scobj.harmony, subset = cell_type == "Epithelial") 
# p_tissue <- VlnPlot(sc.epithelial, c(gene_symbols), group.by = "tissue")
# ggsave(p_tissue, filename = "Vln_tissue_epi(GSE193337).pdf", width=9, height=6, units="in")
# pseudobulk--------------------------------------------------------------------
geom_volcano <- function(dat,pos.num = 10,neg.num = -10, title){
  require(ggplot2)
  require(ggrepel)
  
  dat <- dat %>% 
    tibble::rownames_to_column("SYMBOL") %>% 
    mutate(color = case_when(log2FoldChange >  1 & padj < 0.05 ~ 2,
                             log2FoldChange < -1 & padj < 0.05 ~ 1,
                             TRUE ~ 0)) %>% 
    mutate(label = factor(color, levels = c(0,1,2), labels = c("Non-Sig","Down","Up")),
           color = factor(color))
  dat_up <- dat %>% 
    filter(label == "Up") %>% 
    arrange(desc(log2FoldChange))
  dat_up <- dat_up[1:ifelse(nrow(dat_up) >= 10, 10, nrow(dat_up)),]
  if(is.na(dat_up$SYMBOL[1]) & nrow(dat_up) == 1){
    dat_up[1,] <- c(NA,0,0,0,0,1,1,NA,NA)
    dat_up$log2FoldChange <- dat_up$log2FoldChange %>% as.numeric
    dat_up$padj <- dat_up$padj %>% as.numeric
  }
  dat_down <- dat %>% 
    filter(label == "Down") %>% 
    arrange(log2FoldChange)
  dat_down <- dat_down[1:ifelse(nrow(dat_down) >= 10, 10, nrow(dat_down)),]
  if(is.na(dat_down$SYMBOL[1]) & nrow(dat_down) == 1){
    dat_down[1,] <- c("",0,0,0,0,1,1,NA,NA)
    dat_down$log2FoldChange <- dat_down$log2FoldChange %>% as.numeric
    dat_down$padj <- dat_down$padj %>% as.numeric
  }
  
  ggplot(dat, aes(x = log2FoldChange, y = -log10(padj), col = log2FoldChange, label = SYMBOL)) +
    geom_point(
      # aes(size = !!rlang::sym(abundance))
    )+
    geom_vline(xintercept = c(-1,1), color = "gray80", linetype = 2) +
    geom_hline(yintercept = c(1.30103), color = "gray80", linetype = 2) +
    ylab(expression(-log[10]~(adj.~P~value))) +
    xlab("Log2(Fold Change)") +
    labs(color = "") +
    scale_color_gradient2(high = "red3", mid = "white", low = "blue3",
                          midpoint = 0, na.value = "grey80"
                          #  space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour"
    )+
    scale_size_continuous(range = c(0.1, 4)) +
    geom_text_repel(
      data = dat_up,
      color = "red3",
      size = 5,
      nudge_x = pos.num - as.numeric(dat_up$log2FoldChange),
      segment.size=0.3,
      segment.color="grey",
      direction="y",
      hjust= 0,
      max.overlaps = Inf) +
    geom_text_repel(
      data= dat_down,
      color="blue3",
      size=5,
      nudge_x = neg.num - as.numeric(dat_down$log2FoldChange),
      segment.size = 0.3,
      segment.color = "grey",
      direction="y",
      hjust= 1,
      max.overlaps = Inf) +
    labs(title = title) + 
    # labs(size = expression("Abundance (log2)"),
    #      color = expression("Direction signed"),
    #      title = trait.names[i]) +
    theme_minimal() +
    theme(legend.position = "right",
          legend.title.align = 0, # left align
          legend.title = element_text(margin = margin(t = 15, unit = "pt")) # add more space on top of legend titles
          #legend.spacing.y = unit(1,"cm")
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=20),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          aspect.ratio = 1/1.2, panel.grid.major = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    theme(plot.title = element_text(hjust = 0.5, size=20))
}

# tissue
tmp <- scobj.harmony %>% subset(subset = cell_type == "Epithelial")

counts <- AggregateExpression(tmp, group.by = "orig.ident")
counts <- counts$RNA %>% as.data.frame

group <- ifelse(colnames(counts) %in% c("GSM5793824","GSM5793825","GSM5793826","GSM5793827"), 
                "benign", "tumor")
condition<- factor(group, levels = c("benign","tumor"))
coldata <- data.frame(row.names = colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData = counts %>% as.matrix,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds)  
DEG_t_vs_n <- results(dds, cooksCutoff=FALSE, 
                                     name="condition_tumor_vs_benign", 
                                     independentFiltering = FALSE) %>% 
  as.data.frame %>% 
  arrange(padj)
DEG_t_vs_n <- na.omit(DEG_t_vs_n)

save(counts, file = paste0("pseudobulk_","Epithelial",".Rdata"))
save(DEG_t_vs_n, file = paste0("Epithelial","_DESeq2-DEGs.Rdata"))
write.csv(DEG_t_vs_n, paste0("Epithelial","_DESeq2-DEGs.csv"))





