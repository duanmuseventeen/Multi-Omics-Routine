# Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations
# 10.1038/s41467-018-06318-7
# PMID: 30348985

# GSE115469 five healthy neurologically deceased donors (NDD)

# NDD sc- and sn-RNA-seq raw data and raw count matrices are available at GSE243977, 
# previously published NDD liver data included in this study can be found at GSE115469 
# and GSE185477, and incorporated data from Ramachandran et al.3 can be found at GSE136103.

# Standard scRNA-seq filtering excludes cells with a high ratio of reads from 
# mitochondrial genome transcripts, indicating potential plasma membrane rupture 
# and dissociationbased damage18. This filter is often set at 10%, but hepatocytes
# can have very high mitochondrial content19, thus we chose a threshold of 50% to 
# optimize keeping hepatocytes and removing dead and dying cells. ... As expected, 
# hepatocytes were most susceptible to removal by this filter.

# Further, there are naturally occurring binucleated hepatocytes that are expected 
# to have a high library size, which would make it difficult to distinguish
# between doublets and binucleated cells. A doublet filter set to remove cells with 
# the top 99.9% library size (following standard protocols20), mostly removed 
# hepatocyte (cluster #14) and plasma cell (cluster #7) populations.
# 在技术不成熟的时候，要用一些简单可接受的方法完成实验，而非用严格要求膈应自己

# Importantly, NDD grafts are subjected to the systemic inflammation that 
# accompanies brain death and are thus themselves mildly inflamed23.

# We found six distinct hepatocyte populations that were generally less proliferative cells 
# and showed enriched ALB (Albumin) expression, a hallmark of hepatocytes.

# To study the function of the human clusters, we applied pathway analysis to 
# identify active cellular pathways in each cluster.

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
# GSE115469---------------------------------------------------------------------
Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations

2018
NC

Sonya A MacParland 1 2 3, Jeff C Liu 4, Xue-Zhong Ma 5, Brendan T Innes 4 6, 
Agata M Bartczak 5, Blair K Gage 7, Justin Manuel 5, Nicholas Khuu 8, Juan Echeverri 5, 
Ivan Linares 5, Rahul Gupta 5, Michael L Cheng 9, Lewis Y Liu 10, Damra Camat 5, 
Sai W Chung 10, Rebecca K Seliga 5, Zigong Shao 5, Elizabeth Lee 5, Shinichiro Ogawa 7,
Mina Ogawa 7, Michael D Wilson 6 11, Jason E Fish 9 12, Markus Selzner 5, 
Anand Ghanekar 5, David Grant 5, Paul Greig 5, Gonzalo Sapisochin 5, Nazia Selzner 5, 
Neil Winegarden 8, Oyedele Adeyi 5 9 13, Gordon Keller 7 14 15, Gary D Bader 16 17, 
Ian D McGilvray 18

GSM3178782	Patient 1 Total Liver Homogenate
GSM3178783	Patient 2 Total Liver Homogenate
GSM3178784	Patient 3 Total Liver Homogenate
GSM3178785	Patient 4 Total Liver Homogenate
GSM3178786	Patient 5 Total Liver Homogenate

dat <- data.table::fread("GSE115469/GSE115469_Data.csv.gz")
meta<- data.table::fread("GSE115469/GSE115469_CellClusterType.txt.gz")

all(meta$CellName == colnames(data)[-1])
# [1] TRUE
#### Create Seurat Object---------------------------------------------------------
raw.data <- dat %>% 
  as.data.frame %>% 
  tibble::column_to_rownames("V1") %>% 
  as.matrix %>% 
  Matrix(raw.data, sparse = TRUE) 
GSE115469 <- CreateSeuratObject(counts = raw.data, assay = "RNA", project = "GSE115469")
GSE115469 <- AddMetaData(object = GSE115469, metadata = meta)

scobj <- GSE115469
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
    nFeature_RNA > 500 & 
    nFeature_RNA < 4000 & 
    nCount_RNA < 5000 &
    percent.mt < 10)  

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

scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.3")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "tsne", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.3", reduction = "tsne", label = T)
# Annotation------------------------------------------------------------------
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
  MP <-	c('CD68',	'ITGAM',	'ITGAX',	'HLA+AC0-DRA',	'CSF1R',	'CD14')	# 1
  # plasmacytoid dendritic cell
  pDC <- c('LILRA4',	'CLEC4C',	'GZMB')	# 6 ?9			
  # ILC, innate lymphoid cell
  ILC	<- c('KLRF1',	'KLRC1',	'GZMA',	'GZMB',	'NKG7')	# ?2 4 6 9	
  T_cell <- c('CD3D',	'CD3E',	'CD3G',	'CD8A')	# 2 ?9		
  B_cell <-	c('CD79A',	'CD79B',	'CD19',	'MS4A1') # ?5 7 			
  Plasma <-	c('CD79A',	'IGHA2') # 5 ?7					
  Mast <-	c('KIT',	'TPSAB1',	'TPSB2')				
  Endothelia <-	c('PECAM1',	'CDH5',	'ICAM2',	'KDR',	'ERG') # 3
  Mesenchyme <-	c('PDGFRB',	'ACTA2',	'COL1A1',	'COL1A2',	'COL3A1',	'DES',	'DCN') #12
  Hepatocyte <-	c('ALB',	'TF',	'TTR',	'HNF4A', 'CYP2A6') # 0 10 ?13
  Cholangiocyte <- c('EPCAM',	'KRT19',	'CD24')	# 8 13			
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

VlnPlot(scobj.harmony, features = RBC, group.by = "RNA_snn_res.0.3", layer = "counts", log = TRUE)

scobj.harmony@meta.data$cell_type <- "Unknown"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(0,10)] <- "Hepatocyte"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(1)] <- "Mononuclear phagocytes"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(2)] <- "T cell"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(3)] <- "Endothelia"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(4,6,9)] <- "Lymphoid"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(5)] <- "Plasma" 
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(7)] <- "B cell"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(8,13)] <- "Cholangiocyte"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(11)] <- "RBC"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$RNA_snn_res.0.3 %in% c(12)] <- "Mesenchyme"

p1 <- DimPlot(scobj.harmony, reduction = "umap",
              group.by = "RNA_snn_res.0.3",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj.harmony, reduction = "umap",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj.harmony, reduction = "umap",
              group.by = "orig.ident",
              label = TRUE, pt.size = 0.5) 
p4 <- FeaturePlot(scobj.harmony, reduction = "umap", features = "Targetgene")
p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2)

RidgePlot(object = scobj.harmony, group.by = "cell_type", features = 'Targetgene')
VlnPlot(object = scobj.harmony, group.by = "cell_type", features = 'Targetgene')

save(scobj.harmony, file = "GSE115469-scobj.harmony.annot.Rdata")

