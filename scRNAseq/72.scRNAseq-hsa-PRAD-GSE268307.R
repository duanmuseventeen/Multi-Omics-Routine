# GSM8289374	Localized prostate cancer biopsy, patient 1
# GSM8289375	Localized prostate cancer biopsy, patient 2
# GSM8289376	metastatic hormone naïve prostate cancer biopsy, patient 1
# GSM8289377	metastatic hormone naïve prostate cancer biopsy, patient 2

setwd("GSE268307/")

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
hub <- c("GSM8289374", "GSM8289375", "GSM8289376", "GSM8289377")
names(hub) <- hub

scobj <- hub %>%
  Read10X %>%
  CreateSeuratObject(
    min.cells = 0,
    min.features = 0,
    project = "GSE268307",
    assay = "RNA")

dim(scobj@assays$RNA$counts)
# [1] 36601 52892

scobj <- scobj[!str_detect(rownames(scobj), "^ENSG"),]
scobj <- scobj[!str_detect(rownames(scobj), "^A[CFJLP][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^B[X][0-9]+.[0-9]+$"),]

dim(scobj@assays$RNA$counts)
# [1] 24249 52892

metadata <- scobj@meta.data
metadata <- metadata %>% 
  mutate(
    metastasis = case_when(
    orig.ident %in% c("GSM8289374", "GSM8289375") ~ "Localized",
    orig.ident %in% c("GSM8289376", "GSM8289377") ~ "Metastatic")
    ) %>% as.data.frame
scobj@meta.data <- metadata

table(scobj@meta.data$orig.ident)
# GSM8289374 GSM8289375 GSM8289376 GSM8289377 
# 17190       8054      14805      12843 
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

pdf("run1=QC.pdf")
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA"), ncol = 1) + scale_y_continuous(breaks = c(200, 500, 1000, 4000, 5000,6000,7000,8000,9000,10000))
VlnPlot(scobj, group.by = "orig.ident", features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.rp"), ncol = 1)
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nCount = scobj@meta.data$nCount_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nCount, color = group)) +
  geom_density() +
  scale_x_continuous(limits = c(0,20000)) +
  theme_classic()
dev.off()

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 500 &
    nFeature_RNA < 5000 & 
    nCount_RNA > 500 &
    nCount_RNA < 15000 &
    percent.mt < 20 &
    percent.rp < 40 &
    percent.hb < 1)  

pdf("run1=QC2.pdf")
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA"), ncol = 1) + scale_y_continuous(breaks = c(200, 500, 1000, 4000, 5000,6000,7000,8000,9000,10000))
VlnPlot(scobj, group.by = "orig.ident", features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.rp"), ncol = 1)
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nCount = scobj@meta.data$nCount_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nCount, color = group)) +
  geom_density() +
  scale_x_continuous(limits = c(0,20000)) +
  theme_classic()
dev.off()

table(scobj@meta.data$orig.ident)
# GSM8289374 GSM8289375 GSM8289376 GSM8289377 
# 5189       4859      10057       9843 

dim(scobj)
# [1] 24249 29948
# blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(scobj))
# [1] 239
scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj)
# [1] 24011 29948

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "scobj(GSE268307).Rdata")
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

scobj.harmony <- scobj.harmony.25 %>% 
  RunUMAP(reduction = "harmony", dims = 1:25)
# Visualization-----------------------------------------------------------------
pdf("run1=umap.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "metastasis", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.2", reduction = "umap", label = T)
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
save(scobj.harmony, file = "run1=nf2000_h_pc25 scobj.harmony(GSE268307).Rdata")

scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list
for (i in 1:length(scobj.harmony.split)) {
  # pK Identification (ground-truth) -------------------------------------------
  sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:25)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------
  homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]$RNA_snn_res.0.2)  
  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi <- round(0.1 * nrow(scobj.harmony.split[[i]]@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------
  # scobj.harmony.split[[i]] <- doubletFinder(scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  scobj.harmony.split[[i]] <- doubletFinder(
    scobj.harmony.split[[i]], PCs = 1:25, pN = 0.25, pK = pK, 
    nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
}

save(scobj.harmony.split, file = "run1=scobj.harmony.split(GSE268307).Rdata")

Singlet <- c(
  rownames(scobj.harmony.split[[1]]@meta.data) [scobj.harmony.split[[1]]@meta.data$DF.classifications_0.25 == "Singlet"],
  rownames(scobj.harmony.split[[2]]@meta.data) [scobj.harmony.split[[2]]@meta.data$DF.classifications_0.25 == "Singlet"],
  rownames(scobj.harmony.split[[3]]@meta.data) [scobj.harmony.split[[3]]@meta.data$DF.classifications_0.25 == "Singlet"],
  rownames(scobj.harmony.split[[4]]@meta.data) [scobj.harmony.split[[4]]@meta.data$DF.classifications_0.25 == "Singlet"])

scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)

scobj.harmony@meta.data$dblfinder <- NA
scobj.harmony@meta.data$dblfinder <- "doublet"
scobj.harmony@meta.data$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"

dim(scobj.harmony) # [1] 24011 29948
table(scobj.harmony@meta.data$dblfinder)
# doublet singlet 
# 2542   27406 
# Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1=scobj4run2(GSE268307).Rdata")

# 1000 without regress cell cycle
scobj.harmony.1000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData %>% 
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
pdf("run2=nf1000_h clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony.1000,
     scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
     file = "run2=nf1000_h scobj.harmony.Rdata")

# 1500 without regress cell cycle
scobj.harmony.1500 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData %>% 
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
pdf("run2=nf1500_h clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony.1500,
     scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
     file = "run2=nf1500_h scobj.harmony.Rdata")

# 2000 without regress cell cycle
scobj.harmony.2000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
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
pdf("run2=nf2000_h clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
dev.off()
save(scobj.harmony.2000,
     scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
     file = "run2=nf2000_h scobj.harmony.Rdata")

scobj.harmony <- scobj.harmony.1500.25 %>% 
  RunUMAP(reduction = "harmony", dims = 1:25)
# Annotation====================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.2")

epithelial <- c("EPCAM", "SFN", "KRT5", "KRT8", "KRT14","KLK3") #, "SPRR3"
endothelial <- c("VWF", "PECAM1", "ENG", "CDH5") 
fibroblast <- c("FN1", "DCN", "COL1A1", "COL1A2") 
SMC <- c("TAGLN", "CNN1", "PRKG1", "FOXP2")
pericyte <- c("RGS5", "MCAM", "ACTA2", "MYH11")
T_cell <- c("CD2", "CD3D", "CD3E", "CD3G") 
B_cell <-	c("CD19", "CD79A", "MS4A1", "CD79B") 
plasma <-	c("JCHAIN", "MZB1", "IGHG1", "SDC1") 
# monocyte and macrophage
myeloid <- c("CD68", "CD163", "LYZ", "CD14","FCGR3A","C1QA","C1QB","CST3") 
Neutrophil <- c("CSF3R","CXCL8","G0S2","IFITM2") 
DC <- c('CCR7',	'CLEC9A', 'CD1C',	'IRF7',	'LILRA4')
Mast <-	c('TPSAB1','CPA3','HPGDS','KIT')
Neural <- c('PLP1', 'NRXN1', 'NRXN2', 'NRXN3')
Proliferation <- c("TOP2A","BIRC5","MKI67","DLGAP5")

p <- DotPlot(scobj.harmony, features = c(
  endothelial, epithelial, fibroblast, SMC, pericyte,
  T_cell, B_cell,plasma, myeloid, Neutrophil, Mast , DC, Neural,
  Proliferation
) %>% unique, 
group.by = "RNA_snn_res.0.2") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = "Dotplot marker (GSE268307).pdf", width=9, height=9, units="in")

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(0)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(1)] <- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(2)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(3)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(4)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(5)] <- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(6)] <- "B cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(7)] <- "Mast"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(8)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(9)] <- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(10)]<- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(11)]<- "Neural"

save(scobj.harmony, file = "scobj.harmony(singlet_nf1500_h_noreg_pc25_res0.2)anno.Rdata")

table(scobj.harmony$RNA_snn_res.0.2)
# 0    1    2    3    4    5    6    7    8    9   10   11 
# 5398 5052 3909 2855 2541 2143 1412 1179 1169  985  573  190 table(scobj.harmony@meta.data$cell_type)

table(scobj.harmony$cell_type)
# B cell Endothelial  Epithelial        Mast  Mesenchyme     Myeloid      Neural 
# 1412        4482        5396        1179        8180        1169         190 
# T cell 
# 5398
# Figure2=======================================================================
require(patchwork)
p1 <- DimPlot(scobj.harmony, reduction = "umap", 
              group.by = "cell_type", 
              label = TRUE)
p2 <- DimPlot(scobj.harmony, reduction = "umap", 
              group.by = "RNA_snn_res.1", 
              label = TRUE)
p_out1 <- p1 + p2
ggsave(p_out1, filename = "UMAP(GSE268307).pdf", width=10, height=5, units="in")

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

ggsave(p_fea, filename = "FeaturePlot(GSE268307).pdf", width=8, height=8, units="in")
ggsave(p_vln, filename = "VlnPlot(GSE268307).pdf", width=8, height=8, units="in")

Epi_markers <- FindMarkers(scobj.harmony, group.by = "cell_type", ident.1 = "Epithelial")




