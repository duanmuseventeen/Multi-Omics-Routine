# Targeting Immune-Fibroblast Crosstalk in Myocardial Infarction and Cardiac Fibrosis I----
# GSE217494

# Reference---------------------------------------------------------------------
# https://www.jianshu.com/p/442d890d12ea # about cite-seq
# https://cite-seq.com/ # about cite-seq
# https://www.bilibili.com/video/BV1xA411v7MR/?spm_id_from=333.999.0.0&vd_source=7687d926bc3779b048e880123e66bca4
# https://www.jianshu.com/p/0a28b293d6dd
# An entropy-based metricfor assessing the purity of single cell populations
# https://github.com/PaulingLiu/ROGUE
# https://zhuanlan.zhihu.com/p/667102765
# https://satijalab.org/seurat/articles/integration_introduction
# https://www.jianshu.com/p/7c43dc99c4b1
# https://blog.csdn.net/qq_52813185/article/details/134738735 # graphics
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
# https://cloud.tencent.com/developer/article/1865543
# https://www.jianshu.com/p/d1394ee923a2

rm(list = ls()); gc()
#
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree) # use for determine the optimal resolution
library(ROGUE) # use for determine the optimal resolution
library(harmony)
# Practice Script 1-------------------------------------------------------------
# Practice Script 1: Load data--------------------------------------------------
# Download from # 
hub <- c('sample1' ,
         'sample2' ,
         'sample4' ,
         'sample5' ,
         'sample6' ,
         'sample7' ,
         'sample8' ,
         'sample9' ,
         'sample12',
         'sample13',
         'sample15',
         'sample17',
         'sample27',
         'sample28',
         'sample29',
         'sample30',
         'sample32',
         'sample33',
         'sample34',
         'sample39',
         'sample41',
         'sample42'
)
names(hub) <- hub
sceList.GEX <- lapply(hub, function(x){
  sc10X <- Read10X(x)
  sce <- CreateSeuratObject(
    counts = sc10X$`Gene Expression`,
    project = x,
    min.cells = 0,
    min.features = 0,
    assay = "RNA")
  return(sce)
})

sceList.ADT <- lapply(hub, function(x){
  sc10X <- Read10X(x)
  sce <- CreateSeuratObject(
      counts = sc10X$`Antibody Capture`,
      project = x,
      min.cells = 0,
      min.features = 0,
      assay = "ADT")
  return(sce)
})

# # 10X data contains more than one type and is being returned as a list containing matrices of each type.
# # feature.tsv.gz from cite-seq
# ENSG00000243485 MIR1302-2HG  Gene Expression
# <char>      <char>           <char>
#   1: ENSG00000237613     FAM138A  Gene Expression
# 2: ENSG00000186092       OR4F5  Gene Expression
# 3: ENSG00000238009  AL627309.1  Gene Expression
# 4: ENSG00000239945  AL627309.3  Gene Expression
# 5: ENSG00000239906  AL627309.2  Gene Expression
# ---
#   33812:       HLA-A.B.C   HLA-A.B.C Antibody Capture
# 33813:            CD15        FUT4 Antibody Capture
# 33814:           CD109       CD109 Antibody Capture
# 33815:            CD13       ANPEP Antibody Capture
# 33816:            CD31      PECAM1 Antibody Capture
# 
# # feature.tsv.gz from common scRNAseq
# ENSG00000223972     DDX11L1 Gene Expression
# <char>      <char>          <char>
#   1: ENSG00000227232      WASH7P Gene Expression
# 2: ENSG00000278267   MIR6859-1 Gene Expression
# 3: ENSG00000243485 MIR1302-2HG Gene Expression
# 4: ENSG00000284332   MIR1302-2 Gene Expression
# 5: ENSG00000237613     FAM138A Gene Expression
# ---
#   58228: ENSG00000271254  AC240274.1 Gene Expression
# 58229: ENSG00000275405          U1 Gene Expression
# 58230: ENSG00000275987          U1 Gene Expression
# 58231: ENSG00000277475  AC213203.1 Gene Expression
# 58232: ENSG00000268674     FAM231D Gene Expression

scobj.GEX <- merge(x = sceList.GEX[[1]], y = sceList.GEX[-1], add.cell.ids = hub)
scobj.ADT <- merge(x = sceList.ADT[[1]], y = sceList.ADT[-1], add.cell.ids = hub)

table(scobj.GEX@meta.data$orig.ident)
table(scobj.ADT@meta.data$orig.ident)
# sample1 sample12 sample13 sample15 sample17  sample2 sample27 sample28
# 8484     6383    10053     8410     8034     8865     8546     7325
# sample29 sample30 sample32 sample33 sample34 sample39  sample4 sample41
# 11021     9205    12945     9064    12696     7386     9699     8520
# sample42  sample5  sample6  sample7  sample8  sample9
# 7994     5172    10965     8749     6476    12800
# 
# sample1 sample12 sample13 sample15 sample17  sample2 sample27 sample28
# 8484     6383    10053     8410     8034     8865     8546     7325
# sample29 sample30 sample32 sample33 sample34 sample39  sample4 sample41
# 11021     9205    12945     9064    12696     7386     9699     8520
# sample42  sample5  sample6  sample7  sample8  sample9
# 7994     5172    10965     8749     6476    12800

scobj.GEX[["ADT"]] <- scobj.ADT[["ADT"]]
# Practice Script 2-------------------------------------------------------------
# Practice Script 2: Load data--------------------------------------------------
hub <- c('sample1' ,
         'sample2' ,
         'sample4' ,
         'sample5' ,
         'sample6' ,
         'sample7' ,
         'sample8' ,
         'sample9' ,
         'sample12',
         'sample13',
         'sample15',
         'sample17',
         'sample27',
         'sample28',
         'sample29',
         'sample30',
         'sample32',
         'sample33',
         'sample34',
         'sample39',
         'sample41',
         'sample42'
)
names(hub) <- hub
sc10X <- Read10X(hub)
# sc10X$
#   sc10X$Gene Expression   sc10X$Antibody Capture
# Practice Script 2: Create Seurat Object--------------------------------------------------
sce <- CreateSeuratObject(
  counts = sc10X[["Gene Expression"]],
  min.cells = 0,
  min.features = 0,
  project = "GEX",
  assay = "RNA")

adt <- CreateAssayObject(sc10X[["Antibody Capture"]])
all.equal(colnames(sc10X[["Gene Expression"]]), colnames(sc10X[["Antibody Capture"]]))
# [1] TRUE
sce[["ADT"]] <- adt

table(sce@meta.data$orig.ident)

# sample1 sample12 sample13 sample15 sample17  sample2 sample27 sample28
# 8484     6383    10053     8410     8034     8865     8546     7325
# sample29 sample30 sample32 sample33 sample34 sample39  sample4 sample41
# 11021     9205    12945     9064    12696     7386     9699     8520
# sample42  sample5  sample6  sample7  sample8  sample9
# 7994     5172    10965     8749     6476    12800
# Practice Script 2: Add meta.data----------------------------------------------
meta <- readxl::read_excel("meta.xlsx")
colnames(meta) <- colnames(meta) %>% str_remove("^donor_organism.")
meta <- meta %>% 
  mutate(orig.ident = biomaterial_core.biomaterial_id %>% str_remove("_donor$"))

meta.data <- sce@meta.data %>% 
  tibble::rownames_to_column("barcodes") %>% as.data.frame %>% 
  left_join(meta, by = "orig.ident")
rownames(meta.data) <- meta.data$barcodes

sce@meta.data <- meta.data
# Practice Script 2: QC---------------------------------------------------------
sce[["percent.mt"]] <- PercentageFeatureSet(sce, assay = "RNA", pattern = "^MT-")
sce[["percent.rp"]] <- PercentageFeatureSet(sce, assay = "RNA", pattern = "^RP[SL]")

pdf("P2-QC.pdf")
VlnPlot(sce, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(sce, features = c("nCount_RNA"), ncol = 1)
VlnPlot(sce, features = c("percent.mt"), ncol = 1)
VlnPlot(sce, features = c("percent.rp"), ncol = 1)
FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nFeature = sce@meta.data$nFeature_RNA
) %>% 
  ggplot(aes(x = nFeature)) +
  geom_density() + 
  scale_x_continuous(breaks = c(200,500,1000,3000,4000,5000))
dev.off()

sce <- subset(
  sce, 
  subset = nFeature_RNA > 500 & 
    nFeature_RNA < 6000 & 
    nCount_RNA < 5000 &
    percent.mt < 40 &
    percent.rp < 30)  

pdf("P2-QC-2.pdf")
VlnPlot(sce, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(sce, features = c("nCount_RNA"), ncol = 1)
VlnPlot(sce, features = c("percent.mt"), ncol = 1)
VlnPlot(sce, features = c("percent.rp"), ncol = 1)
dev.off()

# We first perform pre-processing and dimensional reduction on both assays independently. We use standard normalization, but you can also use SCTransform or any alternative method.
DefaultAssay(sce) <- 'RNA'
sce <- sce %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50)

DefaultAssay(sce) <- 'ADT'
sce <- sce %>% 
  NormalizeData(normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% 
  RunPCA(features = rownames(sce), reduction.name = 'apca')

#Examine and visualize PCA results a few different ways
if(TRUE){
  DefaultAssay(sce) <- 'RNA'
  pdf("P2-QC-PCA.pdf")
  print(sce[["pca"]], dims = 1:10, nfeatures = 5)
  VizDimLoadings(sce, dims = 1:2, reduction = "pca")
  DimPlot(sce, reduction = "pca")
  DimHeatmap(sce, dims = 1:6, cells = 500, balanced = TRUE)
  dev.off()
  
  DefaultAssay(sce) <- 'ADT'
  pdf("P2-QC-aPCA.pdf")
  # print(sce[["apca"]], dims = 1:10, nfeatures = 5)
  # VizDimLoadings(sce, dims = 1:2, reduction = "apca")
  # DimPlot(sce, reduction = "apca")
  # DimHeatmap(sce, dims = 1:6, reduction = "apca", cells = 500, balanced = TRUE, ncol = 5)
  
  # sce <- JackStraw(sce, reduction = "apca", dims = 50, num.replicate = 100)
  # # Error in JackStraw(sce, reduction = "apca", dims = 50, num.replicate = 100) :
  # # Only pca for reduction is currently supported
  # sce <- ScoreJackStraw(sce, reduction = "apca", dims = 1:50)
  # JackStrawPlot(sce, reduction = "apca", dims = 1:50)
  ElbowPlot(sce, reduction = "apca", ndims = 50)
  dev.off()
}

DefaultAssay(sce) <- 'RNA'
# Practice Script 2: Not Ingegrate Data-----------------------------------------
if(FALSE){
  # Cluster---------------------------------------------------------------------
  pdf("P2 clustree without integrating data.pdf")
  sce <- FindNeighbors(sce, reduction = "pca", dims = 1:30)
  sce <- FindClusters(sce, resolution = seq(0.1, 1, 0.1))
  clustree(sce, prefix = "RNA_snn_res.")
  dev.off()
  
  # Visualization---------------------------------------------------------------
  #UMAP
  sce <- RunUMAP(sce, reduction = "pca", min_dist = 0.3, dims = 1:30)
  #T-SNE
  sce <- RunTSNE(sce, reduction = "pca", dims = 1:30)
  
  pdf("P2 umap and tsne without integrating data.pdf")
  sce <- RegroupIdents(sce, metadata = "RNA_snn_res.0.3")
  DimPlot(sce, group.by = "orig.ident", reduction = "umap", label = T)
  DimPlot(sce, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T)
  DimPlot(sce, group.by = "orig.ident", reduction = "tsne", label = T)
  DimPlot(sce, group.by = "RNA_snn_res.0.3", reduction = "tsne", label = T)
  dev.off()
}
# Practice Script 2: Ingegrate Data---------------------------------------------
if(FALSE){
  sce.harmony <- RunHarmony(
    sce,
    group.by.vars = "orig.ident",
    reduction.use = "pca",
    reduction.save = "harmony")
  # Practice Script 2: Cluster----------------------------------------------------
  pdf("P2 clustree with integrating data.pdf")
  sce.harmony <- FindNeighbors(sce.harmony, reduction = "harmony", dims = 1:30)
  sce.harmony <- FindClusters(sce.harmony, resolution = seq(0.1, 1, 0.1))
  clustree(sce.harmony, prefix = "RNA_snn_res.")
  dev.off()
  
  # Practice Script 2: Visualization-----------------------------------------------------------------
  #UMAP
  sce.harmony <- RunUMAP(sce.harmony, reduction = "harmony", min_dist = 0.3, dims = 1:30)
  # head(sce.all@reductions$umap@cell.embeddings) 
  #note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
  # DimPlot(sce.harmony, group.by = "RNA_snn_res.0.4", reduction = "umap", label = T)
  #T-SNE
  sce.harmony <- RunTSNE(sce.harmony, reduction = "harmony", dims = 1:30)
  # head(sce.all@reductions$tsne@cell.embeddings)
  # DimPlot(sce.harmony, group.by = "RNA_snn_res.0.4", reduction = "tsne", label = T) 
  
  pdf("P2 umap and tsne with integrating data.pdf")
  sce.harmony <- RegroupIdents(sce.harmony, metadata = "RNA_snn_res.0.3")
  DimPlot(sce.harmony, group.by = "orig.ident", reduction = "umap", label = T)
  DimPlot(sce.harmony, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T)
  DimPlot(sce, group.by = "orig.ident", reduction = "tsne", label = T)
  DimPlot(sce, group.by = "RNA_snn_res.0.3", reduction = "tsne", label = T)
  
  sce.harmony <- RegroupIdents(sce.harmony, metadata = "RNA_snn_res.0.4")
  DimPlot(sce.harmony, group.by = "orig.ident", reduction = "umap", label = T)
  DimPlot(sce.harmony, group.by = "RNA_snn_res.0.4", reduction = "umap", label = T)
  DimPlot(sce, group.by = "orig.ident", reduction = "tsne", label = T)
  DimPlot(sce, group.by = "RNA_snn_res.0.4", reduction = "tsne", label = T)
  dev.off()
}
# Practice Script 2: WNN--------------------------------------------------------
if(TRUE){
  # Identify multimodal neighbors. These will be stored in the neighbors slot, 
  # and can be accessed using sce[['weighted.nn']]
  # The WNN graph can be accessed at sce[["wknn"]], 
  # and the SNN graph used for clustering at bm[["wsnn"]]
  # Cell-specific modality weights can be accessed at sce$RNA.weight
  
  sce <- FindMultiModalNeighbors(
    sce, 
    reduction.list = list("pca", "apca"), 
    dims.list = list(1:30, 1:20), 
    modality.weight.name = "RNA.weight"
  )
  # Warning message:
  # In FindMultiModalNeighbors(sce, reduction.list = list("pca", "apca"),  :
  # The number of provided modality.weight.name is not equal to the number of modalities. RNA.weight ADT.weight are used to store the modality weights
  
  sce <- FindClusters(sce, graph.name = "wsnn", algorithm = 3, resolution = seq(0.1,1,0.1), verbose = FALSE)
  
  sce <- RunUMAP(sce, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  
  pdf("clustree WNN.pdf")
  clustree(sce, prefix = "wsnn_res.")
  dev.off()
  
  pdf("RNA vs WNN.pdf")
  p1 <- DimPlot(sce, reduction = 'wnn.umap', group.by = 'wsnn_res.0.6', label = TRUE, 
                repel = TRUE, label.size = 2.5) + NoLegend()
  p2 <- DimPlot(sce, reduction = 'umap', group.by = 'RNA_snn_res.0.6', label = TRUE, 
                repel = TRUE, label.size = 2.5) + NoLegend()
  p1 + p2
  dev.off()
  
  pdf("RNA vs ADT.pdf")
  sce <- RunUMAP(sce, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                 reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
  sce <- RunUMAP(sce, reduction = 'apca', dims = 1:20, assay = 'ADT', 
                 reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
  p3 <- DimPlot(sce, reduction = 'rna.umap', group.by = 'wsnn_res.0.6', label = TRUE, 
                repel = TRUE, label.size = 2.5) + NoLegend()
  p4 <- DimPlot(sce, reduction = 'adt.umap', group.by = 'wsnn_res.0.6', label = TRUE, 
                repel = TRUE, label.size = 2.5) + NoLegend()
  p3 + p4
  dev.off()
}
# Practice 2: Cell Cycle (Only observe)--------------------------------------------------------
# https://satijalab.org/seurat/articles/cell_cycle_vignette
# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(
  file = "nestorawa_forcellcycle_expressionMatrix.txt",
  header = TRUE, as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sce <- CellCycleScoring(sce, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Visualize the distribution of cell cycle markers across
pdf("cell_cycle.pdf")
RidgePlot(sce, group.by = "Phase", features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
sce <- RunPCA(sce, features = c(s.genes, g2m.genes), reduction.name = "pca.cell.cycle")

pdf("cell_cycle2.pdf")
DimPlot(sce, reduction = "pca.cell.cycle", group.by = "Phase", label = TRUE, pt.size = 0.5) 
DimPlot(sce, reduction = "umap", group.by = "Phase", label = TRUE, pt.size = 0.5) 
dev.off()

# if batch effect by cell cycle
# sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(sce))
# or 
# sce$CC.Difference <- sce$S.Score - sce$G2M.Score
# sce <- ScaleData(sce, vars.to.regress = "CC.Difference", features = rownames(sce))
# Practice 2: Annotation--------------------------------------------------------------------
# wsnn_res.0.6
sce <- RegroupIdents(sce, metadata = "wsnn_res.0.6")
sce.copy <- sce

ACM <- c('TNNT2', 'TTN', 'MYH6', 'MYH7') # 30
VCM <- c('MYH7', 'MYL2', 'FHL2') # 30
EC <- c('VWF', 'CDH5', 'THBD', 'TEK', 'ENG', 'PECAM1') # 1 12 13 15 19 2 21 26 3 5 7 8
FB <- c('FN1', 'VIM', 'COL1A1', 'COL1A2') # 0 10 14 16 18 23 ?26 27 ?28 31 4 ?9
M  <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74') # 17 20 22 23 24 6
Myeloid	<- c('CD163', 'CD14', 'C1QA') # 17 20 22 23 24 ?27 6
P <- c('HCN1', 'HCN4', 'CACNA1D') # pacemaker cells (P cells) 
Peri <- c('RGS5', 'ABCC9', 'KCNJ8') # 14 16 4 9
SMC <- c('MYH11', 'ACTA2','CNN1', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') # 14 16 4 9
Meso <- c('MSLN', 'WT1', 'BNC1') # 
Neural <- c('PLP1', 'NRXN1', 'NRXN3') # 28
Adi <- c('GPAM', 'FASN', 'LEP') # 
Lymph <- c('CD3E', 'IL7R', 'CD40LG') # 11 18 ?25
Mast <- c('KIT', 'CPA3') # 29
NKT <- c('SLAMF1', 'SLAMF6', 'TGFBR', 'VA24', 'JA18') # 
RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ')

if(TRUE){
  # pdf("Vln.pdf")
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = ACM, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = VCM, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = EC, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = FB, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = M, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = Myeloid, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = P, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = Peri, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = SMC, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = Meso, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = Neural, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = Adi, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = Lymph, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = Mast, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = NKT, layer = "counts", log = TRUE)
  VlnPlot(sce, group.by = "wsnn_res.0.6", features = RBC, layer = "counts", log = TRUE)
  # dev.off()
}

sce@meta.data$cell_type <- "Unknown" # 32 (86 cells)
sce@meta.data$cell_type[sce@meta.data$wsnn_res.0.6 == 30] <- "CM"
sce@meta.data$cell_type[sce@meta.data$wsnn_res.0.6 %in% c(1, 2, 3, 5, 7, 8, 12, 13, 15, 19, 21, 26)] <- "EC"
sce@meta.data$cell_type[sce@meta.data$wsnn_res.0.6 %in% c(0, 4, 10, 14, 16, 18, 23, 27, 31)] <- "FB"
sce@meta.data$cell_type[sce@meta.data$wsnn_res.0.6 %in% c(6, 17, 20, 22, 23, 24, 27)] <- "Myeloid"
sce@meta.data$cell_type[sce@meta.data$wsnn_res.0.6 %in% c(4, 9, 14, 16)] <- "SMC" # "Peri"
sce@meta.data$cell_type[sce@meta.data$wsnn_res.0.6 == 28] <- "Neural"
sce@meta.data$cell_type[sce@meta.data$wsnn_res.0.6 %in% c(11, 18, 25)] <- "Lymphoid"
sce@meta.data$cell_type[sce@meta.data$wsnn_res.0.6 == 29]<- "Mast"

p1 <- DimPlot(sce, reduction = "umap",
              group.by = "wsnn_res.0.6",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(sce, reduction = "umap",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(sce, reduction = "umap",
              group.by = "orig.ident",
              label = TRUE, pt.size = 0.5) 
p1 + p2 + p3
# ggsave(..., width = 18, height = 6)

save(sce, file = "sce.annot.Rdata")









