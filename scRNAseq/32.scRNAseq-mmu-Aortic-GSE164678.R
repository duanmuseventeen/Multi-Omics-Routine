# 	Single-cell RNAseq analysis of CaCl2-induced mouse abdominal aortic aneurysm (AAA) samples----
# GSE164678

# Summary	We performed single-cell RNA sequencing (scRNA-seq) on infrarenal abdominal aortas from C57BL/6J mice after perivascular CaCl2 treatment. Infrarenal abdominal aortas were collected four days after AAA induction and processed for sequencing. These data provide high-resolution insight into the complexity and heterogeneity of mouse AAA.
# Overall design	AAA was induced in 12-week-old male C57BL/6J mice. Infrarenal abdominal aortas were perivascularly treated with 0.5 M CaCl2 for 10 minutes followed by PBS for 5 minutes. Mice in the sham group received 0.5 M NaCl treatment for 10 minutes followed by PBS for 5 minutes. Aortas were collected four days after AAA induction. Single-cell suspensions from four aortas were pooled as one sample.

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
library(AUCell)
# Load data---------------------------------------------------------------------
hub <- c("abdominal_aortic_aneurysm(AAA)", "sham")
names(hub) <- hub

scobj.list <- lapply(hub, function(x){
  x %>% 
    Read10X %>% 
    CreateSeuratObject(
      counts = ., 
      project = x,
      min.cells = 0,
      min.features = 0,
      assay = "RNA") %>% 
    return
})

scobj.list

scobj <- merge(x = scobj.list[[1]], y = scobj.list[-1])

scobj

table(scobj@meta.data$orig.ident)
# abdominal_aortic_aneurysm(AAA)                           sham 
# 1922                           3322 
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^Mt")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^Rp[sl]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^Hb[^(p)]")

VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nFeature = scobj@meta.data$nFeature_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nFeature, color = group)) +
  geom_density() +
  scale_x_continuous(breaks = c(200,300,400,500,750,1000,3000,4000,5000))

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 100 &
    nFeature_RNA < 7500 & 
    nCount_RNA < 20000 &
    percent.mt < 10 &
    percent.rp < 20 &
    percent.hb < 1)  

table(scobj@meta.data$orig.ident)
# abdominal_aortic_aneurysm(AAA)                           sham 
# 1168                           2064 

VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)

scobj <- scobj %>% 
  SCTransform %>% # vst.flavor = 'v2', verbose = FALSE
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

clustree(scobj.harmony, prefix = "SCT_snn_res.")
# Visualization-----------------------------------------------------------------
#UMAP
scobj.harmony <- RunUMAP(scobj.harmony, reduction = "harmony", min_dist = 0.3, dims = 1:30)
#T-SNE
scobj.harmony <- RunTSNE(scobj.harmony, reduction = "harmony", dims = 1:30)

DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "SCT_snn_res.0.6", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "SCT_snn_res.0.6", reduction = "tsne", label = T)
# Filter doublets---------------------------------------------------------------
if(TRUE){
  # scDblFinder
  scobj.harmony.sce <- as.SingleCellExperiment(scobj.harmony, assay = "SCT")
  scobj.harmony.sce <- scDblFinder(
    scobj.harmony.sce, 
    samples = scobj.harmony.sce$orig.ident,
    dbr.sd=1)
  
  scobj.harmony <- as.Seurat(scobj.harmony.sce)
  scobj.h.sc <- subset(scobj.harmony, scDblFinder.class == "singlet")
}
# Annotation--------------------------------------------------------------------
scobj.h.sc <- RegroupIdents(scobj.h.sc, metadata = "SCT_snn_res.0.6")

SMC <- c('MYH11', 'ACTA2','CNN1', 'SMTN', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') 
# Peri <- c('RGS5', 'ABCC9', 'KCNJ8')
FB <- c('LAMA2', 'TSHZ2', "FN1", "COL1A1") 
EC <- c('PECAM1',"VWF") 
RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ') # 
Meso <- c('MSLN', 'WT1', 'BNC1') 
Adi <- c('GPAM', 'FASN', 'LEP') 
Neural <- c('PLP1', 'NRXN1', 'NRXN3') 
Lymph <- c('CD3E', 'IL7R', 'CD40LG')
Myeloid <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74', 'CD14', 'C1QA')
Mast <- c("GATA2","HPGDS","HPGD","KIT","IL1RL1","LAPTM4A","MLPH","VWA5A","RGS13","LTC4S") 
NK <- c("XCL2","XCL1","KLRD1","GNLY","KLRC1","CD7","CTSW","PRF1","CMC1")
plasma <- c("IGHA1","MZB1","SSR4","PRDX4","FKBP11","SEC11C","QPCT","LMF1")
B <- c("CD19","MS4A1","LY9","BANK1","ARHGAP24","CD79B","LTB","LINC00926","STAG3")
Mast <- c("GATA2","HPGDS","HPGD","KIT","IL1RL1","LAPTM4A","MLPH","VWA5A","RGS13","LTC4S") # 14

if(TRUE){
  # Single-cell suspensions were prepared from sham and CaCl2-treated infrarenal 
  # aortas 4 days after AAA induction. Verhoeff Van Gieson elastic stain and 
  # immunostaining of contractile SMC marker MYH11 (SMC-specific myosin heavy chain) 
  # on aortic cross sections showed significant elastin degradation and SMC 
  # depletion at this time point (Supplemental Figure I). A total of 3896 cells 
  # were recovered after application of quality control filters, including 2537 
  # cells from the sham group, and 1359 cells from the AAA group. Unsupervised 
  # Seurat-based clustering identified 12 distinct cell populations (Figure 1A). 
  # We assigned putative biological identities to each population based on gene 
  # expression patterns of established canonical markers of SMCs (Myh11, Acta2), 
  # endothelial cells (Cdh5, Pecam1), fibroblasts (Col1a1), macrophages (Lyz2), 
  # neutrophils (S100a8 and S100a9), dendritic cells (Klrd1 and Flt3), T & natural 
  # killer cells (Cd3g, Gzma), and B cells (Cd79a, Ms4a1)4, 5, 16 (Figure 1B and 
  # Supplemental Figure II).
  
  SMCs <- c('Myh11', 'Acta2','Tagln','Aqp1') # 0 4
  fibroblasts <- c(
    'Dcn', 'Sfrp4', 'Lum', 'Col1a1',
    'Col1a1', 'Pclaf') # ?0 1 ?4 ?5
  EC <- c('Cdh5', 'Cldn5', 'Vwf', 'Pecam1') # 8 10
  macrophages <- c(
    'Pf4','Lyz2', 'Mrc1',
    'Lyz2', 'Il1b', 'Mrc1',
    'Stmn1') # 2 3 6 9
  neutrophils <- c('S100a8', 'S100a9','Cxcr4','Cxcl2','Il1b','Mmp9','Cd14') # 7
  DC <- c('Klrd1', 'Flt3', 'H2-Ab1') # 11
  T_and_NK <- c('Cd3g', 'Gzma') # 11
  B_cell <- c('Cd79a', 'Ms4a1', 'Cd19') # ?6
  cell_cycle <- c("Mki67", "Ccna2", "Ccnb1", "Pcna", "Top2a", "Mcm6")
}

VlnPlot(scobj.h.sc, features = fibroblasts, group.by = "SCT_snn_res.0.6", layer = "counts", log = TRUE)

scobj.h.sc@meta.data$cell_type <- "Unknown"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.6 %in% c(0, 4)] <- "SMC"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.6 %in% c(1, 5)] <- "FB"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.6 %in% c(2, 3, 6, 9)] <- "macrophages"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.6 %in% c(9)] <- "cell_cycle"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.6 %in% c(7)] <- "neutrophils"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.6 %in% c(8, 10)] <- "EC"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.6 %in% c(11)] <- "DC & T & NK"

p1 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "SCT_snn_res.0.6",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "orig.ident",
              label = TRUE, pt.size = 0.5) 
p1 + p2 + p3
# ggsave(..., width = 18, height = 6)

FeaturePlot(object = scobj.h.sc, slot = "counts", features = 'Targetgene')
RidgePlot(object = scobj.h.sc, features = 'Targetgene')
DotPlot(scobj.h.sc, group.by = "disease", features = "Targetgene")
DotPlot(scobj.h.sc, group.by = "cell_type", features = "Targetgene")

save(scobj.h.sc, file = "GSE164678-scobj.annot.Rdata")
