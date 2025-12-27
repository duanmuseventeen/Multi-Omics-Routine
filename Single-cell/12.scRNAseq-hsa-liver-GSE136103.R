# Resolving the fibrotic niche of human liver cirrhosis at single-cell level
# PMID:31597160
# 10.1038/s41586-019-1631-3

# GSE136103

rm(list = ls());gc()
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree) # use for determine the optimal resolution
library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
# Load data---------------------------------------------------------------------
hub <- c(
  "GSM4041151_healthy1_cd45-A",
  "GSM4041154_healthy2_cd45",
  # "GSM4041157_healthy3_cd45-B",
  # "GSM4041152_healthy1_cd45-B",
  "GSM4041156_healthy3_cd45-A",
  "GSM4041159_healthy4_cd45-")
names(hub) <- hub

sc10X <- Read10X(hub)

GSE136103 <- CreateSeuratObject(
  counts = sc10X,
  min.cells = 0,
  min.features = 0,
  project = "GSE136103",
  assay = "RNA")

scobj <- GSE136103

table(scobj@meta.data$orig.ident)
# GSM4041151 GSM4041152 GSM4041154 GSM4041156 GSM4041157 GSM4041159
# 1039        478       4688       3230       1156       3306

table(scobj@meta.data$orig.ident)
GSM4041151 GSM4041154 GSM4041156 GSM4041159
1039       4688       3230       3306
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
    nFeature_RNA > 400 & 
    nFeature_RNA < 4000 & 
    nCount_RNA < 25000 &
    percent.mt < 10)  

pdf("QC.pdf")
VlnPlot(scobj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
dev.off()

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

pdf("clustree.pdf")
clustree(scobj.harmony, prefix = "RNA_snn_res.")
dev.off()
# Visualization-----------------------------------------------------------------
#UMAP
scobj.harmony <- RunUMAP(scobj.harmony, reduction = "harmony", min_dist = 0.3, dims = 1:30)
#T-SNE
scobj.harmony <- RunTSNE(scobj.harmony, reduction = "harmony", dims = 1:30)

pdf("umap and tsne.pdf")
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.4")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "tsne", label = T)
dev.off()
# Annotation--------------------------------------------------------------------
# LSEC
# The liver vascular endothelium, made up of LSECs and the endothelium of blood
# vessels, provides a dynamic barrier between the blood and the liver
# microenvironment. LSECs are defined as scavenger endothelia that use clathrin-mediated 
# endocytosis for the clearance of macromolecules from the blood36.

# Cholangiocytes. 
# Cholangiocytes are epithelial cells which line bile ducts that comprise 3–5% of 
# the cells in the human liver and generate 40% of the total bile volume40. 
# Intrahepatic cholangiocytes originate from hepatoblasts during embryonic development.
# Based on the duct diameter size, cholangiocytes are categorized as small and large, 
# each with different secretory functions and sensitivity to injury20. The circumference 
# of small bile ducts is formed by 4–5 cholangiocytes and by 10–12 cholangiocytes in
# larger ducts41.

# Hepatic stellate cells. 
# Hepatic stellate cells (HSCs) are found in the subendothelial space of the liver 
# sinusoid, known as the space of Disse46. HSCs are the main storage site for Vitamin A 
# and are major contributors to tissue fibrosis. Upon activation, human HSCs express 
# α-smooth muscle actin (ACTA2) and begin to lay down extracellular matrix, which is 
# composed of collagen (e.g., COL1A1, COL1A2).

# We observed three endothelial cell populations, which were less proliferative than immune cells
# 因此，细胞周期的不同是细胞亚群的不同的一种体现，对细胞周期进行调整并不十分合理

# GSE243977
if(TRUE){
  liverMarker <- list(
    `B cell` = c("CD79A","CD79B"),
    `Naive B cell` = c("IGHM"),
    `Plama B cell` = c("IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGKC",
                       "IGLC2","IGLC3","IGLL5"),
    cellCycle = c("BIRC5","CDC20","CDK1","CDK4","HMGB2","MKI67","TOP2A"),
    ApoLipoChol = c("APOA1", "APOA2", "APOC1", "APOC3"),
    Cholangio = c("AGXT","AMBP","ANXA4","CLDN1","DEFB1","EEF1A1","EPCAM","FXYD2",
                  "GNB2L1","HP","MT2A","NEAT1","RPL3"),
    KeratinsChol = c("KRT18","KRT19","KRT7","KRT8"),
    Mucus Secreting = c("LGALS2","MUC1","MUC3A","MUC5B","PIGR","SCGB3A1""SLPI","SPINK1","TFF3"),
    Arterial = c("CD34","PLVAP","PODXL"),
    cvEndo = c("ACKR1","RSPO3","WNT2"),
    cvLSEC = c("CLEC1B","CLEC4G","CLEC4M","CRHBP","FCN2","FCN3","STAB1","STAB2"),
    ppLSEC = c("ADRIF","AQP1","CLEC14A","ID1","IGFBP7","MGP","S100A6","SPARCL1","TM4SF1"),
    Endothelial = c("C7","DNASE1L3","ENG","GSN","INMT","LIFR","PECAM1","PLAC8","RAMP3","TIMP3","VWF"),
    `C-Hepato` = c("ADH1A","ADH1B","ADH1C","ADH4","CES1","CYP2E1","CYP3A4","CYP3A5","DCXR","GSTA1","OAT"),
    `P-Hepato` = c("AGT","ALB","ALDOB","APOA1","CYP1A2","CYP2A6","CYP2A7","CYP2D6","FGA","FGB","GLS2","HAL","HAMP","SDS"),
    Hepatocyte = c("ASGR1","SERPINA1"),
    `T cell` = c("CD3D","CD3E","TRAC","TRBC2","CD4","CD40LG","IL7R","LTB","RCAN3","S100A4","SPOCK2","CD8A","CD8B"),
    cNK = c("FCGR3A","GNLY","GZMB","TBX21","XCL2"),
    gdT = c("CD7","CTSW","GATA3","TRDC","TRGC1","TRGC2"),
    lrNK = c("CD69","CMC1","EOMES","GZMK","KLRC1","KLRF1","XCL1"),
    ActivatedMac = c("AREG","CCL3","CD83","CXCL2","CXCL3","IL1B","NAMPT","PLAUR","SRGN","THBS1"),
    Kupffer = c("C1QA","C1QB","C1QC","CD163","CD5L","CTSB","FTL","HMOX1","MARCO","SEPP1","SLC40A1","VCAM1"),
    Macrophage = c("APOE","CD14","CD54","CD68","FTH1","IL18","NINJ1","PLAC8","PLTP","VSIG4"),
    `LAM-like` = c("ACP5","CSTB","FABP5","LGMN","PLD3","PSAP"),
    `MHCII/cDC` = c("CD74","HLA-DPA1","HLA-DQB1","HLA-DRA","HLA-DRB1"),
    `Monocyte-derived` = c("FCN1","LYZ","MNDA","S100A12","S100A4","S100A6","S100A8","S100A9","VCAN"),
    Myeloid = c("CST3"),
    Neutrophil = c("BASP1","CSF3R","CXCL8","CXCR1","CXCR2","FCGR3B","FPR1","G0S2","IFITM2","S100A11"),
    ResidentMac = c("FCER1G","FOLR2","LYVE1","MS4A7","TIMD4","TIMP1"),
    SynapseMac = c("AIF1","COTL1","IFITM3","LST1"),
    RBC = c("HBA1","HBA2","HBB","HBD"),
    `AP1+ Stellate` = c("CCL2","FBLN5","FOSB","JUNB","PDGFRA","SOCS3","SOD2"),
    `Quiescent Stellate` = c("COLEC11","HGF","PTH1R","RBP1","RELN","VIPR1"),
    `Stellate-Fiber` = c("ACTA2","COL1A1","COL1A2","MYL9","TAGLN"),
    VSMC = c("CD9","CRYAB","PMEPA1","RGCC")
  )
}
# GSE136103
if(TRUE){
  # MP, mononuclear phagocytes
  MP <-	c('CD68',	'ITGAM',	'ITGAX',	'HLA+AC0-DRA',	'CSF1R',	'CD14')	# 1 3
  # plasmacytoid dendritic cell
  pDC <- c('LILRA4',	'CLEC4C',	'GZMB')				
  # ILC, innate lymphoid cell
  ILC	<- c('KLRF1',	'KLRC1',	'GZMA',	'GZMB',	'NKG7') # 4 5 14 ?15		
  T_cell <- c('CD3D',	'CD3E',	'CD3G',	'CD8A') # 0 2			
  B_cell <-	c('CD79A',	'CD79B',	'CD19',	'MS4A1') # 12 
  Plasma <-	c('CD79A',	'IGHA2') # 13	
  Mast <-	c('KIT',	'TPSAB1',	'TPSB2')	# ?0			
  Endothelia <-	c('PECAM1',	'CDH5',	'ICAM2',	'KDR',	'ERG')	# 6 8 9 11 17	
  Mesenchyme <-	c('PDGFRB',	'ACTA2',	'COL1A1',	'COL1A2',	'COL3A1',	'DES',	'DCN') # 10 16
  Hepatocyte <-	c('ALB',	'TF',	'TTR',	'HNF4A',	'CYP2A6')	# 7	
  Cholangiocyte <- c('EPCAM',	'KRT19',	'CD24')	# 7			
  RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ')
  CellCycle	<- c('MKI67', 'CCNA2',	'CCNB2',	'STMN1')			
  inflammatory <- c("S100A8", "S100A9", "S100A6", "VCAN", "LYZ")
  non_inflammatory <- c("CD5L", "MARCO", "VCAM1") 
}
# GSE185477
if(TRUE){
  RBC  <- c("HBB","HBA1","HBA")
  CD3T <- c("TRAC","TRBC1","TRBC2","TRDC","TRGC1","TRGC2","CD8A","CD8B","CD3D","CD3G","IL32")
  BCR  <- c("IGKC","IGHE","IGHM","IGLC1","IGLC3","IGLC2","JCHAIN","IGHA1")
  MATB <- c("CD22","CD37","CD79B","FCRL1","LTB","DERL3","IGHG4")
  CLOT <- c("F10","F5","F9","F7","FGB")
  BILE <- c("SLC10A1","SLC01BA","SLC01B3","SLC22A1","SLC22A7")
  HEP  <- c("CYP2E1","CYP3A7","CYP2A7","CYP2A6","CYP1A2","CYP3A4","GLUL","DCXR","FTL","GPX2","GSTA1","FABP1","HAL","AGT","ALDOB","SDS")
  LSEC <- c("FCN2","CLEC1B","CLEC4G","PVALB","S100A13","GJA5","SPARCL1","CLEC14A","PLVAP","EGR3")
  NKT  <- c("CSTW","IL7R","GZMB","GZMH","TBX21","HOPX","PRF1","S100B","TRDC","TRGC1","TRGC2","IL2RB","KLRB1","NCR1","NKG7","NCAM1","XCL2","XCL1","CD160","KLRC1")
  MAC  <- c("VCAN","S100A8","MNDA","LYZ","FCN1","CXCL8","VCAM1","TTYH3","TIMD4","SLC40A1","RAB31","MARCO","HMOX1","C1QC")
  CHOL <- c("PROM1","SOX9","KRT7","KRT19","CFTR","EPCAM","CLDN4","CLD7","ANXA4","TACSTD2")
  STEL <- c("ACTA2","COL1A1","RBP1","TAGLN","ADAMTSL2","GEM","LOXL1","LUM")
  ENDO <- c("PECAM1","TAGLN","VWF","FLT1","MMRN1","RSPO3","LYPD2","LTC4S","TSHZ2","IL1R1")
}
# GSE115469
if(TRUE){
  hept <- c("CYP3A7","CYP2A7","CYP2A6",
            "SCD","HMGCS1","ACSS2","TM7SF2","SEC16B","SLBP","RND3","PCK1","BCHE","G6PC","GHR","ALDH6A1","RPP25L","HSD11B1","HAMP","GHR","HPR","GSTA2","AKR1C1","MASP2")  
  lsec <- c("SPARCL1","CLEC14A","S100A13","CLEC1B","FCN2",
            "CCL14","MGP","TM4SF1")
  ec <- c("RAMP3","INMT","DNASE1L3","LIFR")
  chol <- c("KRT7","KRT19","SOX9","EPCAM")
  stel <- c("ACTA2","COL1A1","RBP1")
  macr <- c("S100A8","LYZ","S100A9","HLA-DPB1","CD5L","MARCO","VSIG4")
  t_cell <- c("CD2","CD3D","TRAC","GZMK","GNLY","PTGDS","GZMB","TRDC","STMN1","HMGB2","TYMS")
  b_cell <- c("MS4A1","LTB","CD37","CD79B")
  plasma <- c("IGLC2","IGHG1","IGKC")
  nk <- c("CD7","KLRB1","NKG7")
  erthyroid <- c("HBB","CA1","ALAS2")
}

VlnPlot(scobj.harmony, features = RBC, group.by = "RNA_snn_res.0.4", layer = "counts", log = TRUE)

scobj.harmony@meta.data$cell_type <- "Unknown"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.4 %in% c(0,2)] <- "T cell"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.4 %in% c(1,3)] <- "Mononuclear phagocytes"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.4 %in% c(4,5,14,15)] <- "Lymphoid"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.4 %in% c(6,8,9,11,17)] <- "Endothelia"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.4 %in% c(7)] <- "Hepatocyte and Cholangiocyte"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.4 %in% c(10,16)] <- "Mesenchyme" 
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.4 %in% c(12)] <- "B cell" 
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.4 %in% c(13)] <- "Plasma" 


p1 <- DimPlot(scobj.harmony, reduction = "umap",
              group.by = "RNA_snn_res.0.4",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj.harmony, reduction = "umap",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj.harmony, reduction = "umap",
              group.by = "orig.ident",
              label = TRUE, pt.size = 0.5) 
p1 + p2 + p3

FeaturePlot(scobj.harmony, reduction = "umap", features = "Targetgene")
VlnPlot(scobj.harmony, group.by = "cell_type", features = "Targetgene")
RidgePlot(scobj.harmony, group.by = "cell_type", features = "Targetgene")

save(scobj.harmony, file = "GSE136103-scobj.harmony.annot.Rdata")






