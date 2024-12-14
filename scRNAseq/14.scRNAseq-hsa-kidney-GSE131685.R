# Single-Cell RNA Sequencing in Multiple Pathologic Types of Renal Cell Carcinoma Revealed Novel Potential Tumor-Specific Markers
# Single-cell RNA sequencing of human kidney

# GSE131685

# normal adjacetn tissue

# Fresh human kidney samples (Supplementary Table S1) were collected at the First 
# Affiliated Hospital of Guangxi Medical University and Affiliated Tumour Hospital of 
# Guangxi Medical University. Two out of the three samples were obtained from patients 
# undergoing radical nephrectomy, and the remaining sample was from a patient undergoing 
# radical nephroureterectomy. Normal kidney tissues were obtained at least 2 cm away 
# from tumour tissue.

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
# Load data---------------------------------------------------------------------
GSM4145204	Patient 1 Kidney Homogenate
GSM4145205	Patient 2 Kidney Homogenate
GSM4145206	Patient 3 Kidney Homogenate

hub <- dir()
names(hub) <- hub
scobj <- hub %>% 
  Read10X %>% 
  CreateSeuratObject(
    project = "GSE131685",
    min.cells = 0,
    min.features = 0,
    assay = "RNA")

table (scobj@meta.data$orig.ident)
# GSM4145204 GSM4145205 GSM4145206 
# 8164       6499      10741 
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nFeature = scobj@meta.data$nFeature_RNA,
  group = scobj@meta.data$orig.ident
) %>% 
  ggplot(aes(x = nFeature, color = group)) +
  geom_density() + 
  scale_x_continuous(breaks = c(200,300,400,500,1000,3000,4000,5000))

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 300 & 
    nFeature_RNA < 2000 & 
    nCount_RNA < 6000 &
    percent.mt < 30 &
    percent.rp < 40 &
    percent.hb < 1)  

VlnPlot(scobj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)

scobj <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50)
# Integration-------------------------------------------------------------------
scobj.harmony <- RunHarmony(
  scobj,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  reduction.save = "harmony")
# Cluster-----------------------------------------------------------------------
scobj.harmony <- FindNeighbors(scobj.harmony, reduction = "harmony", dims = 1:30)
scobj.harmony <- FindClusters(scobj.harmony, resolution = seq(0.1, 1, 0.1))

clustree(scobj.harmony, prefix = "RNA_snn_res.")
# Visualization-----------------------------------------------------------------
#UMAP
scobj.harmony <- RunUMAP(scobj.harmony, reduction = "harmony", min_dist = 0.3, dims = 1:30)
#T-SNE
scobj.harmony <- RunTSNE(scobj.harmony, reduction = "harmony", dims = 1:30)

DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T)
# DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.6", reduction = "umap", label = T)
# DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.7", reduction = "umap", label = T)
# Cell cycle--------------------------------------------------------------------
exp.mat <- read.table(file = "/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.harmony <- CellCycleScoring(scobj.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(scobj.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

scobj.harmony <- RunPCA(scobj.harmony, features = c(s.genes, g2m.genes))
DimPlot(scobj.harmony)
# Filter doublets---------------------------------------------------------------
# Filter doublets & Investigate batch effect by Shannon Index
if(FALSE){
  scobj.harmony.sce <- as.SingleCellExperiment(scobj.harmony)
  scobj.harmony.sce <- scDblFinder(scobj.harmony.sce, dbr.sd=1)
  
  # observed. We quantified the batch effect using Shannon entropy
  # as described in the CellMixS R package47. Mean entropy of the transcript data
  # was 0.91 and 0.92 for ADT data (Supplementary Fig. 4c,d). We observed that cells
  # from a single cluster (indicated by the black arrow) had lower entropy values and
  # contributed most to the batch effect.
  scobj.harmony.sce <- CellMixS::entropy(scobj.harmony.sce, group = "orig.ident", k = 20)
  
  scobj.harmony <- as.Seurat(scobj.harmony.sce)
  scobj.h.sc <- subset(scobj.harmony, scDblFinder.class == "singlet")
  
  # scobj.h.sc@meta.data %>% 
  #   mutate(mean = mean(entropy)) %>% 
  #   head
  # 
  # scobj.h.sc@meta.data %>% 
  #   group_by(RNA_snn_res.0.3) %>% 
  #   mutate(mean = mean(entropy)) %>% 
  #   distinct(mean)
}
# Annotation--------------------------------------------------------------------
# curated on the basis of 35657798
if(TRUE){
  # proximal tubular cells
  PT <- c("MIOX", "ALDOB", "FABP1", "PCK1", "ANPEP") # 0 2 3 ?4
  # tubular progenitor cells (PG)
  PG <- c("CD24", "PROM1", "MMP7") # 8
  SMC <- c("ACTA2", "RGS5", "MYH11", "TAGLN") # 
  Macrophage <- c("CD68", "CD163", "LYZ") # 6
  Monocyte <- c("CD14", "LYZ", "S100A12", "S100A9", "S100A8") # 6
  DC <- c("FCER1A", "CD1E", "CD1C", "HLA-DMA", "HLA-DMB") # 6 ?9
  NK <- c("KLRD1", "KLRC1", "GZMB", "PRF1") # 5
  Fibroblast <- c("SFRP2", "SPARC", "MMP2", "COL3A1", "COL1A1", "COL1A2", 
                  "EMILIN1", "PDGFRB")
  EC <- c("PECAM1", "PLVAP", "CDH5", "KDR")
  CD8T <- c("CD3D", "CD3E", "CD8A") # 
  CD4T <- c("CD3E", "CD3D", "IL7R") # 
  T_cell <- c("CD3E", "CD3D", "IL7R", "GZMK")# 5
  B <- c("CD79A", "CD79B", "MS4A1", "CD19") # 9
  Plasma <- c("IGKC") # 9
  Mast <- c("TPSAB1", "TPSB2", "KIT")
  # CellCycle	<- c('MKI67', 'CCNA2', 'CCNB2', 'PCNA', 'STMN1')
  RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ')
}
# 37468583
if(TRUE){
  ####epithelial cells----
  renal_corpuscle <- c('PTPRQ', 'WT1', 'NTNG1', 'NPHS1', 'NPHS2', 'CLIC5', 'PODXL')
  glomerulus <- c('CLDN1', 'VCAM1', 'CFH', 'RBFOX1', 'ALDH1A2') 
  proximal_tubules <- c(
    'LRP2', 'CUBN', 'SLC13A1', 'SLC5A12', 'SLC13A3', 'SLC22A6', 'PRODH2', 'SLC5A2', 
    'SLC22A8', 'SLC34A1', 'SLC22A7', 'MOGAT1', 'SLC5A11', 'SLC22A24', 'SLC7A13', 
    'SLC5A8', 'ABCC3', 'SATB2'
  ) # 0 2 3 4 ?1 
  PT <- c("MIOX", "ALDOB", "FABP1", "PCK1", "ANPEP")
  intermediate_tubules <- c(
    'CRYAB', 'TACSTD2', 'SLC44A5', 'KLRG2', 'COL26A1', 'BOC', 'VCAM1', 'SLC39A8', 
    'AQP1', 'LRRC4C', 'LRP2', 'UNC5D', 'SATB2', 'JAG1', 'ADGRL3', 'ID1', 'CLDN1', 
    'AKR1B1', 'CLDN4', 'BCL6', 'SH3GL3', 'SLC14A2', 'SMOC2', 'BCAS1', 'CLCNKA', 
    'CLDN10', 'PROX1'
  ) #8
  Distal_tubules <- c(
    'CASR', 'SLC12A1', 'UMOD',
    'NELL1', 'ESRRB', 'EGF', 'CLDN14', 'PROX1', 'MFSD4A', 'KCTD16', 'RAP1GAP', 'ANK2', 
    'CYFIP2', 'PPM1E', 'GP2', 'ENOX1', 'TMEM207', 'TMEM52B', 'CLDN16', 'WNK1M',
    'NOS1', 'ROBO2', 'CALCR', 'PPFIA2', 'PAPPA2', 'SLC12A3', 'CNNM2', 'FGF13', 'KLHL3', 
    'LHX1', 'TRPM6', 'TRPM7', 'ADAMTS17', 'ITPKB', 'ZNF385D', 'HS6ST2',
    'TRPV5', 'SLC8A1', 'SCN2A', 'HSD11B2', 'CALB1'
  ) %>% unique #7
  Collecting_tubules <- c(
    'SLC8A1', 'SCN2A', 'HSD11B2', 'CALB1', 'KITLG', 'PCDH7', 'RALYL', 'TOX', 'SGPP1', 
    'SCNN1G', 'SCNN1B', 'KCNIP1', 'GATA3', 'AQP2', 'AQP3', 'FXYD4', 'SOX5', 'PDE10A', 
    'SLC25A29', 'ST6GAL1', 'PAPPA', 'SYK', 'FAM81A', 'PROM1', 'KCNK13', 'FXYD4', 'SOX5',
    'PHACTR1', 'PCDH7', 'SLC14A2', 'HS3ST5', 'TACSTD2', 'TP63', 'GPX2', 'FXYD3', 'KRT5',
    'ATP6V0D2', 'ATP6V1C2', 'TMEM213', 'CLNK', 'SLC4A1', 'SLC26A7', 'HS6ST3', 'NXPH2', 
    'LEF1', 'ADGRF5', 'SLC8A1', 'SCN2A', 'CALB1', 'KIT', 'AQP6', 'STAP1', 'FAM184B', 
    'CALCA', 'SLC4A9', 'SLC35F3', 'SLC26A4', 'INSRR', 'TLDC2'
  ) %>% unique
  #### Endothelial Cell----
  #### Vascular Smooth Muscle Cell / Pericyte(stroma cells)----
  #### Fibroblast(stroma cells)----
  #### Immune Cells----
  B Cell <- BANK1, BLK, MS4A1, BACH2
  Plasma Cell <- IGKC, TENT5C, MZB1, FCRL5, CD38, JCHAIN
  T Cell <- CD96, CD247, THEMIS, BCL11B, CAMK4, IL7R
  Natural Killer T Cell <- CD96, CD247, RUNX3, GNLY, NKG7, CCL5, KLRF1, CCL4, GZMA
  Mast Cell <- MS4A2, CPA3, KIT
  M2 Macrophage <- F13A1, MRC1, CD163, STAB1, SLC1A3, CD14, FOLR2
  Monocyte-derived Cell <- MSR1, ITGAX, HLA-DQA1, HLA_DRB1, CSF2RA, CD14, TRPM2
  Classical Dendritic Cell <- ITGAX, HLA-DQA1, HLA-DRA, CSF2RA, CIITA, WDFY4, FLT3, ZNF366, CADM1, ZBTB46, CLEC9A
  Plasmacytoid Dendritic Cell <- IRF8, CUX2, P2RY14, IL3RA, CLEC4C
  Non-Classical Monocyte <- CTSS, IRAK3, TCF7L2, TNFRSF1B, FCN1, HLA-DRA, FCGR3A
  Neutrophil <- S100A9, S100A8, IFITM2, FCGR3B, CD1C
}

scobj.h.sc <- RegroupIdents(scobj.h.sc, "RNA_snn_res.0.3")
VlnPlot(scobj.h.sc, features = renal_corpuscle, group.by = "RNA_snn_res.0.3", layer = "counts", log = TRUE)

scobj.h.sc@meta.data$cell_type <- "Unknown"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.3 %in% c(0,2,3,4)] <- "Proximal tubules"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.3 %in% c(1)] <- "Proximal tubules"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.3 %in% c(5)] <- "T cell" 
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.3 %in% c(6)] <- "Myeloid"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.3 %in% c(7)] <- "Distal tubules"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.3 %in% c(8,10)] <- "Intermediate tubules & Collecting tubules"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.3 %in% c(9)] <- "B cell and plasma"

p1 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "RNA_snn_res.0.3",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "orig.ident",
              label = TRUE, pt.size = 0.5) 
p4 <- FeaturePlot(scobj.h.sc, reduction = "UMAP", features = "Targetgene")
p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2)

RidgePlot(object = scobj.h.sc, group.by = "cell_type", features = "Targetgene")
VlnPlot(object = scobj.h.sc, group.by = "cell_type", features = "Targetgene")

save(scobj.h.sc, file = "GSE131685-scobj.h.sc.annot.Rdata")