Single-cell RNA sequencing based on Seq-Well S^3 protocol of 
  6 prostate biopsies from 3 different patients, 
  4 radical prostatectomies with tumor-only samples from 4 patients, and 
  4 radical prostatectomies with matched normal samples from 4 patients. 
Organoids derived from radical prostatectomy paired tumor and normal tissue samples.

# Given limited gene symbols, this dataset is not recommended to be used.

setwd("GSE176031_Seq-well/")

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

setwd("rawdat/")
hub        <- dir()
names(hub) <- hub
mtx        <- lapply(hub, function(x){
    tmp <- data.table::fread(x)
    colnames(tmp)[-1] <- paste0(str_remove(x, "_.*"), "-", colnames(tmp)[-1])
    return(tmp)
  }) 
names(mtx) <- names(mtx) %>% str_remove_all("_.*")
lapply(mtx, dim)
lapply(mtx, function(x){ sum(duplicated(x$GENE)) }) %>% unlist %>% unique
commongenes<- Reduce(intersect, lapply(mtx, function(x){x$GENE}))
length(commongenes) # [1] 10799
mtx        <- lapply(mtx, function(x){
  x %>% filter(GENE %in% commongenes) %>% as.data.frame %>% tibble::column_to_rownames("GENE")
  })
lapply(mtx, dim)
lapply(mtx, function(x){ rownames(x) }) %>% Reduce(intersect, .) %>% length
mtx.integ  <- bind_cols(mtx)
dim(mtx.integ) # [1] 10799 34845

setwd("GSE176031_Seq-well/")
metadata <- data.frame(
  orig.ident = str_remove_all(colnames(mtx.integ), "-.*"),
  row.names = colnames(mtx.integ)) %>% 
  mutate(geo_accession = orig.ident) %>% 
  left_join(readxl::read_excel("GSE176031_series_matrix.xlsx"), by = "geo_accession")
# Create Seurat object----------------------------------------------------------
scobj <- CreateSeuratObject(
  counts = mtx.integ %>% as.matrix %>% Matrix(sparse = TRUE), 
  meta.data = metadata, 
  project = "GSE176031")

scobj
# An object of class Seurat 
# 10799 features across 34845 samples within 1 assay 
# Active assay: RNA (10799 features, 0 variable features)
# 1 layer present: counts

scobj <- scobj[!str_detect(rownames(scobj), "^A[CFJLP][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^B[X][0-9]+.[0-9]+$"),]
scobj <- scobj[!str_detect(rownames(scobj), "^ENSG"),]

dim(scobj)
# [1] 10799 34845

table(scobj@meta.data$orig.ident)
# GSM5353214 GSM5353215 GSM5353216 GSM5353217 GSM5353218 GSM5353219 GSM5353220 GSM5353221 
# 766        808        952        827        946        978        917       1092 
# GSM5353222 GSM5353223 GSM5353224 GSM5353225 GSM5353226 GSM5353227 GSM5353228 GSM5353229 
# 832        647        854       1184       1179       1029        632        714 
# GSM5353230 GSM5353231 GSM5353232 GSM5353233 GSM5353234 GSM5353235 GSM5353236 GSM5353237 
# 1053       1142       1295       1151        742        816        935       1104 
# GSM5353238 GSM5353239 GSM5353240 GSM5353241 GSM5353242 GSM5353243 GSM5353244 GSM5353245 
# 1191       1099       1305       1039        956       1127       1152       1148 
# GSM5353246 GSM5353247 GSM5353248 
# 1110       1068       1055 
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

pdf("run1=QC.pdf")
# VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA"), ncol = 1) + NoLegend()
VlnPlot(scobj, group.by = "orig.ident", features = c("nCount_RNA"), ncol = 1) + NoLegend()
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20)) + NoLegend()
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.rp"), ncol = 1) + scale_y_continuous(breaks = c(10,20,30,40)) + NoLegend()
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.hb"), ncol = 1) + NoLegend()
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 500 &
    nFeature_RNA < 3000 &
    nCount_RNA > 500 &
    nCount_RNA < 5000 &
    percent.mt < 20 &
    percent.rp < 10 &
    percent.hb < 1)  

pdf("run1=QC2.pdf")
VlnPlot(scobj, group.by = "orig.ident", features = c("nFeature_RNA"), ncol = 1) + NoLegend()
VlnPlot(scobj, group.by = "orig.ident", features = c("nCount_RNA"), ncol = 1) + NoLegend()
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20)) + NoLegend()
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.rp"), ncol = 1) + NoLegend()
VlnPlot(scobj, group.by = "orig.ident", features = c("percent.hb"), ncol = 1) + NoLegend()
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

table(scobj$orig.ident)
# GSM5353214 GSM5353215 GSM5353216 GSM5353217 GSM5353218 GSM5353219 GSM5353220 GSM5353221 
# 225        364        612        398        492        513        553        605 
# GSM5353222 GSM5353223 GSM5353224 GSM5353225 GSM5353226 GSM5353227 GSM5353228 GSM5353229 
# 510        318        365        850        797        642        127         42 
# GSM5353230 GSM5353231 GSM5353232 GSM5353233 GSM5353234 GSM5353235 GSM5353236 GSM5353237 
# 679        672        239        319        328        399        528        624 
# GSM5353238 GSM5353239 GSM5353240 GSM5353241 GSM5353242 GSM5353243 GSM5353244 GSM5353245 
# 500        245        380        387        251        566        552         68 
# GSM5353246 GSM5353247 GSM5353248 
# 473        467        373

dim(scobj)
# [1] 10799 15463
# blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("blacklist.xlsx")
sum(blacklist$Symbol %in% rownames(scobj))
# [1] 190

scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
dim(scobj)
# [1] 10610 15463

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "run1=scobj(GSE176031).Rdata")
# Nomalization & harmony & cluster----------------------------------------------
scobj.harmony <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "patient"),
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
pdf("run1=nf2000_h_pc25.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T) + NoLegend()
DimPlot(scobj.harmony, group.by = "patient", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "disease", reduction = "umap", label = T)
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
save(scobj.harmony, file = "run1=scobj.harmony(GSE176031).Rdata")

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

save(scobj.harmony.split, file = "run1=scobj.harmony.split(GSE176031).Rdata")

Singlet <- c()
for (i in 1:length(scobj.harmony.split)) {
  Singlet <- c(Singlet, 
               rownames(scobj.harmony.split[[i]]@meta.data) [scobj.harmony.split[[i]]@meta.data$DF.classifications_0.25 == "Singlet"])
}
for (i in 1:length(scobj.harmony.split)) {
  all(scobj.harmony.split[[i]]@meta.data[rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet,   54] == "Singlet") %>% print
  all(scobj.harmony.split[[i]]@meta.data[!(rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet),54] != "Singlet") %>% print
}

scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)

# scobj.h.dblfinder <- subset(scobj.harmony, subset = id %in% Singlet)
scobj.harmony@meta.data$dblfinder <- NA
scobj.harmony@meta.data$dblfinder <- "doublet"
scobj.harmony@meta.data$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"

dim(scobj.harmony) # [1] 10610 15463
table(scobj.harmony@meta.data$dblfinder)
# doublet singlet 
# 649   14814
# Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1=scobj4run2(GSE176031).Rdata")

# 1000----
scobj.harmony.1000 <- scobj %>% 
  subset(subset = dblfinder == "singlet") %>%
  NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData(vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), 
            features = rownames(scobj)) %>% 
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = c("orig.ident", "patient"),
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
    group.by.vars = c("orig.ident", "patient"),
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
    group.by.vars = c("orig.ident", "patient"),
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


scobj.harmony <- scobj.harmony.1000.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Annotation====================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.5")

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
group.by = "RNA_snn_res.0.5") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = "Dotplot marker (GSE176031).pdf", width=9, height=8, units="in")

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(0)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(1)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(2)] <- "T cell"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(3)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(4)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(5)] <- "Epithelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(6)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(7)] <- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(8)] <- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(9)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(10)]<- "Doublets"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(11)]<- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(12)]<- "Proliferation"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(13)]<- "Endothelial" # Mesenchyme
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(14)]<- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.5 %in% c(15)]<- "Myeloid"

# C10 <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.5", ident.1 = 10)
# C11 <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.5", ident.1 = 11)
# C14 <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.5", ident.1 = 14)
# C15 <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.5", ident.1 = 15)

# save(scobj.harmony, file = "scobj.harmony(singlet_nf1000_h_regccmt_pc30_res0.5)anno.Rdata")
load("scobj.harmony(singlet_nf1000_h_regccmt_pc30_res0.5)anno.Rdata")

table(scobj.harmony$RNA_snn_res.0.5)
# 0    1    2    3    4    5    6    7    8    9   10   11 
# 3316 3004 2066 1683 1278 1154  899  492  469  279  239  223 
# 12   13   14   15 
# 115   98   88   60 
table(scobj.harmony$cell_type)
# Doublets   Endothelial    Epithelial    Mesenchyme       Myeloid Proliferation 
# 239           997          8752          1049          2022           115 
# T cell       Unknown 
# 2066           223 
# Figure2=======================================================================
require(patchwork)
p1 <- DimPlot(scobj.harmony, reduction = "umap", 
              group.by = "cell_type", 
              label = TRUE)
p2 <- DimPlot(scobj.harmony, reduction = "umap", 
              group.by = "RNA_snn_res.1", 
              label = TRUE)
p_out1 <- p1 + p2
ggsave(p_out1, filename = "UMAP(GSE176031).pdf", width=10, height=5, units="in")

p_fea <- lapply(gene_symbols, 
                function(x){FeaturePlot(scobj.harmony, features = x) + 
                    scale_colour_gradientn(
                      colours = colorRampPalette(
                        c('gray90','#FFCA62','#FFB336','#FF9700','#FF5A00','#F24410',
                          '#E52C22','#DD1D23','#C20030','#930039','#8C003A',
                          '#6F003D','#56033F'))(1000))}) %>% 
  patchwork::wrap_plots(ncol = 2, nrow = 2)
p_vln <- VlnPlot(scobj.harmony, features = gene_symbols, 
                 group.by = "cell_type")

ggsave(p_fea, filename = "FeaturePlot(GSE176031).pdf", width=8, height=8, units="in")
ggsave(p_vln, filename = "VlnPlot(GSE176031).pdf", width=8, height=8, units="in")
