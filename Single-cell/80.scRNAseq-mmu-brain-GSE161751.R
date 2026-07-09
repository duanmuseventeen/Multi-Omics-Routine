rm(list = ls());gc()

setwd("GSE161751_RAW/")
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
# library(infercnv) # CNV
# library(copykat) # CNV

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

colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

hub <- dir()
names(hub) <- hub

gse161751 <- hub %>%
  Read10X %>%
  CreateSeuratObject(
    min.cells = 0,
    min.features = 0,
    project = "gse161751",
    assay = "RNA")

gse161751$tissue   <- "PVN"
gse161751$group    <- NA
gse161751$group[gse161751$orig.ident == "GSM4914022"] <- "Control"
gse161751$group[gse161751$orig.ident == "GSM4914023"] <- "Stress"
# I - run1|QC gse161751 -------------------------------------------------------------
setwd("/")
gse161751[["percent.mt"]] <- PercentageFeatureSet(gse161751, pattern = "^mt-")
gse161751[["percent.rp"]] <- PercentageFeatureSet(gse161751, pattern = "^Rp[sl]")
gse161751[["percent.hb"]] <- PercentageFeatureSet(gse161751, pattern = "^Hb[^(p)]")

pdf("run1 gse161751_QC前.pdf")
myqc4seurat(seurat.obj = gse161751)
dev.off()

gse161751.qc <- subset(
  gse161751, 
  subset = 
    nFeature_RNA > 1000 &
    nFeature_RNA < 6000 & 
    nCount_RNA > 500 &
    nCount_RNA < 20000 &
    percent.mt < 10 &
    percent.rp < 15 &
    percent.hb < 1)  

pdf("run1 gse161751_QC后.pdf")
myqc4seurat(seurat.obj = gse161751.qc)
dev.off()

table(gse161751.qc$orig.ident)
# GSM4914022 GSM4914023 
# 3024       2903

scobj <- gse161751.qc
# I - run1|raw data stat--------------------------------------------------------
gtf_data <- rtracklayer::import("refdata-gex-mm10-2020-A_genes.gtf") %>% as.data.frame
keep <- gtf_data$gene_name[gtf_data$gene_type != "lncRNA"] %>% unique

scobj  <- subset(scobj , features = keep)

dim(scobj)
# [1] 20713  5927
# I - run1|blacklist------------------------------------------------------------
# blacklist <- readxl::read_excel("blacklist.xlsx")
# blacklist <- blacklist %>% filter(complete.cases(Mouse))

# sum(blacklist$Mouse %in% rownames(scobj))
# # [1] 23
# scobj <- scobj[!(rownames(scobj) %in% blacklist$Symbol),]
# dim(scobj)
# # [1] 19006 45820
# 
# scobj <- JoinLayers(scobj, assay = "RNA")
scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)

save(scobj, file = "run1 scobj.Rdata")
# I - run1|Nomalization & harmony & cluster----------------------------------------------
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

pdf("run1 nf2000_h_clustree 20 25 30.pdf")
ElbowPlot(scobj.harmony, reduction = "pca", ndims = 50)
clustree(scobj.harmony.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.30, prefix = "RNA_snn_res.")
dev.off()

scobj.harmony <- scobj.harmony.20 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20)
# I - run1|Visualization-----------------------------------------------------------------
pdf("run1 nf2000_h_pc20 umap.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.6", reduction = "umap", label = T)
FeaturePlot(scobj.harmony, features = "percent.mt")
dev.off()
# I - run1|Cell Cycle--------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# 将人源基因转换为匹配的鼠源基因
library(biomaRt)
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

conversion_table <- getLDS(
  attributes = c("hgnc_symbol"),
  filters = "hgnc_symbol",
  values = c(s.genes, g2m.genes),
  mart = human_mart,
  attributesL = c("mgi_symbol"),
  martL = mouse_mart,
  uniqueRows = TRUE 
)

s.genes.m   <- conversion_table$MGI.symbol[conversion_table$HGNC.symbol %in% s.genes]
g2m.genes.m <- conversion_table$MGI.symbol[conversion_table$HGNC.symbol %in% g2m.genes]

scobj.harmony <- CellCycleScoring(scobj.harmony, s.features = s.genes.m, g2m.features = g2m.genes.m, set.ident = TRUE)

pdf("run1 nf2000_h_pc20 umap cellcycle.pdf")
RidgePlot(scobj.harmony, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
DimPlot(scobj.harmony, reduction = "umap")
dev.off()
# I - run1|doublets----------------------------------------------------------------------
save(scobj.harmony, file = "run1 nf2000_h_pc20 scobj.harmony.Rdata")

scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list

nPC = 20
for (i in 1:length(scobj.harmony.split)) {
  # pK Identification (ground-truth) -------------------------------------------
  sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:nPC)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
  ## Homotypic Doublet Proportion Estimate -------------------------------------
  homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]$RNA_snn_res.0.5)  
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

save(scobj.harmony.split, file = "run1 nf2000_h_pc20 scobj.harmony.split.Rdata")

Singlet <- c()
for (i in 1:length(scobj.harmony.split)) {
  Singlet <- c(Singlet, 
               rownames(scobj.harmony.split[[i]]@meta.data) [scobj.harmony.split[[i]]@meta.data$DF.classifications_0.25 == "Singlet"])
}
finalcol <- ncol(scobj.harmony.split[[1]]@meta.data)
for (i in 1:length(scobj.harmony.split)) {
  all(scobj.harmony.split[[i]]@meta.data[rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet,   25] == "Singlet") %>% print
  all(scobj.harmony.split[[i]]@meta.data[!(rownames(scobj.harmony.split[[i]]@meta.data) %in% Singlet),25] != "Singlet") %>% print
}

scobj.harmony$id <- rownames(scobj.harmony@meta.data)

# scobj.h.dblfinder <- subset(scobj.harmony, subset = id %in% Singlet)
scobj.harmony$dblfinder <- NA
scobj.harmony$dblfinder <- "doublet"
scobj.harmony$dblfinder[scobj.harmony$id %in% Singlet] <- "singlet"

dim(scobj.harmony)
# [1] 20713  5927
table(scobj.harmony$dblfinder)
# doublet singlet 
# 271    5656 

pdf("run1 nf2000_h_pc20 doublets.pdf")
DimPlot(scobj.harmony, group.by = "dblfinder", reduction = "umap", label = T)
dev.off()
# I - run2|Re-run========================================================================
scobj[['S.Score']] <- scobj.harmony$S.Score
scobj[['G2M.Score']] <- scobj.harmony$G2M.Score
scobj[['Phase']] <- scobj.harmony$Phase
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="run1 scobj4run2.Rdata")
# I - run2|noreg================================================================
# I - run2|1000 without regress cell cycle----
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
pdf("run2 nf1000_h_noreg pc20-30 clustree.pdf")
ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
dev.off()
# I - run2|1500 without regress cell cycle----
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
pdf("run2 nf1500_h_noreg pc20-30 clustree.pdf")
ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
dev.off()
# I - run2|2000 without regress cell cycle----
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
pdf("run2 nf2000_h_noreg pc20-30 clustree.pdf")
ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
dev.off()
# I - run2|Run Umap----
scobj.harmony <- scobj.harmony.1000.20 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20)
# I - run2|Visualization========================================================
pdf("run2 nf1500_h_pc30_regccmt_res1.0 umap.pdf")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T) + NoLegend() 
DimPlot(scobj.harmony, group.by = "RNA_snn_res.1", reduction = "umap", label = T)
dev.off()
# I - run2|Annotation====================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.1")

# Reference: 
# 1. Nat Neurosci, 2017 PMID: 27991900
# 2. Nat Neurosci, 2017 PMID: 28166221
# 3. NC, 2019 PMID: 31444346

neuron = c("Tubb3","Avp","Oxt","Rgs16","Tac2","Ghrh","Slc18a2")
VLMC   = c("Col1a1","Col3a1","Lum")
NG2    = c("Cspg4")                   # OPC
macro  = c("Aif1")                    # PVMmacro
endo   = c("Slco1c1")
mural  = c("Mustn1")
cellcyle  = c("Top2a", "Mki67", "Birc5", "Pcna")
astrocyte = c("Gfap")
parstuber = c("Cyp2f2", "Tshb", "Timeless", "Cck")
tanycyte  = c("Adm","Rax","Crym")
ependymocyte = c("Ccdc153")
oligodendrocyte = c("Mag")

GABAergic = c("Gad1","Gad2","Slc32a1")
glutamate = "Slc17a6" # glutamatergic
dopamine  = c("Th","Slc6a3")

p_dotplot <- DotPlot(scobj.harmony, features = c(
  neuron, VLMC, NG2, macro, endo, mural, cellcyle, astrocyte, parstuber,
  ependymocyte, tanycyte, oligodendrocyte,
  GABAergic, glutamate, dopamine
) %>% unique, 
group.by = "RNA_snn_res.1") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(p_dotplot, filename = "DotPlot4Anno(regccmt res_1).pdf", width=10, height=12, units="in")

Clist <- list()
n <- 1
for (i in 0:20) {
  Clist[[n]] <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.1", ident.1 = as.character(i))
  write.csv(Clist[[n]], paste0("run2 nf1000_h_noreg_pc20_res1_C",n-1,".csv"))
  
  n <- n + 1
}

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(0)] <- "Oligodendrocyte"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(1)] <- "GABAergic Neuron"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(2)] <- "Astrocyte"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(3)] <- "glutamatergic Neuron"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(4)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(5)] <- "Oligodendrocyte"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(6)] <- "Ependymocyte"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(7)] <- "GABAergic Neuron"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(8)] <- "GABAergic Neuron"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(9)] <- "GABAergic Neuron"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(10)]<- "Pars Tuber"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(11)]<- "Macrophage"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(12)]<- "NG2"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(13)]<- "NG2"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(14)]<- "Tanycyte"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(15)]<- "Macrophage"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(16)]<- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(17)]<- "GABAergic Neuron"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(18)]<- "glutamatergic Neuron"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(19)]<- "Mural"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.1 %in% c(20)]<- "Cycling"

table(scobj.harmony$RNA_snn_res.1)
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 842 580 525 459 389 357 354 294 271 241 232 227 172 154 153 131  94  61  61  44  15 

table(scobj.harmony$cell_type)
# Astrocyte              Cycling          Endothelial         Ependymocyte 
# 525                   15                  483                  354 
# GABAergic Neuron glutamatergic Neuron           Macrophage                Mural 
# 1447                  520                  358                   44 
# NG2      Oligodendrocyte           Pars Tuber             Tanycyte 
# 326                 1199                  232                  153

# save(scobj.harmony, file = "scobj.harmony(singlet_nf1500_h_regccmt_pc30_res0.7)anno.Rdata")
load("scobj.harmony(singlet_nf1000_h_noreg_pc20_res1)anno.Rdata")
# I - Visualization=============================================================
fig <- scobj.harmony
fig$sample <- fig$orig.ident
fig$RNA_snn_res.1 <- factor(fig$RNA_snn_res.1, levels = c(0:20))

fig$group[fig$orig.ident == "GSM4914022"] = "Control"
fig$group[fig$orig.ident == "GSM4914023"] = "Stress"
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

focus = readxl::read_excel("泛素化酶.xlsx")

library(biomaRt)
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

conversion_table <- getLDS(
  attributes = c("hgnc_symbol"),
  filters = "hgnc_symbol",
  values = focus$`gene symbol`,
  mart = human_mart,
  attributesL = c("mgi_symbol"),
  martL = mouse_mart,
  uniqueRows = TRUE 
)

setwd("泛素化相关基因")
genes = conversion_table$MGI.symbol[conversion_table$MGI.symbol %in% rownames(fig)]
genes = genes[!duplicated(genes)]
for (gene in genes) {
  pdf(paste0(gene,".pdf"))
  VlnPlot(fig, features = gene, group.by = "cell_type") %>% print
  dev.off()
}

Clist = list()
celltypes = unique(fig$cell_type)  
for (celltype in celltypes) {
  Clist[[celltype]] <- FindMarkers(scobj.harmony, group.by = "cell_type", ident.1 = celltype)
  write.csv(Clist[[celltype]], paste0(celltype,".csv"))
}
  
  
  
  
  
  
  
  
  
