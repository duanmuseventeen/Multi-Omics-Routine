# Single-nucleus Ribonucleic Acid-sequencing and Spatial Transcriptomics Reveal the Cardioprotection of Shexiang Baoxin Pill (MUSKARDIA) in Mice with Myocardial Ischemia-Reperfusion Injury
# PMID: 37229263
# 10.3389/fphar.2023.1173649

# GSE227088

GSM7091453	sham, replicate 1, scRNAseq
GSM7091454	sham, replicate 2, scRNAseq
GSM7091455	sham, replicate 3, scRNAseq
GSM7091456	I/R, replicate 1, scRNAseq
GSM7091457	I/R, replicate 2, scRNAseq
GSM7091458	I/R, replicate 3, scRNAseq
GSM7091459	MUSKARDIA, replicate 1, scRNAseq
GSM7091460	MUSKARDIA, replicate 2, scRNAseq
GSM7091461	MUSKARDIA, replicate 3, scRNAseq

rm(list = ls())
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree) # use for determine the optimal resolution
library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
# linux-------------------------------------------------------------------------
for i in $(seq 3 8);
do
mkdir GSM709145${i}
mv GSM709145${i}*_barcodes.tsv.gz GSM709145${i}/barcodes.tsv.gz
mv GSM709145${i}*_features.tsv.gz GSM709145${i}/features.tsv.gz
mv GSM709145${i}*_matrix.csv.gz   GSM709145${i}/matrix.csv.gz
done
# Load data---------------------------------------------------------------------
GSM7091453 <- data.table::fread("GSM7091453/matrix.csv.gz")   
GSM7091454 <- data.table::fread("GSM7091454/matrix.csv.gz")   
GSM7091455 <- data.table::fread("GSM7091455/matrix.csv.gz")   
GSM7091456 <- data.table::fread("GSM7091456/matrix.csv.gz")   
GSM7091457 <- data.table::fread("GSM7091457/matrix.csv.gz")   
GSM7091458 <- data.table::fread("GSM7091458/matrix.csv.gz") 

colnames(GSM7091453) <- paste0(str_replace(colnames(GSM7091453), "\\.", "-"), "-GSM7091453")
colnames(GSM7091454) <- paste0(str_replace(colnames(GSM7091454), "\\.", "-"), "-GSM7091454")
colnames(GSM7091455) <- paste0(str_replace(colnames(GSM7091455), "\\.", "-"), "-GSM7091455")
colnames(GSM7091456) <- paste0(str_replace(colnames(GSM7091456), "\\.", "-"), "-GSM7091456")
colnames(GSM7091457) <- paste0(str_replace(colnames(GSM7091457), "\\.", "-"), "-GSM7091457")
colnames(GSM7091458) <- paste0(str_replace(colnames(GSM7091458), "\\.", "-"), "-GSM7091458")

colnames(GSM7091453)[1] <- "GENE"
colnames(GSM7091454)[1] <- "GENE"
colnames(GSM7091455)[1] <- "GENE"
colnames(GSM7091456)[1] <- "GENE"
colnames(GSM7091457)[1] <- "GENE"
colnames(GSM7091458)[1] <- "GENE"

scobj <- GSM7091453 %>% 
  full_join(GSM7091454, by = "GENE") %>% 
  full_join(GSM7091455, by = "GENE") %>% 
  full_join(GSM7091456, by = "GENE") %>% 
  full_join(GSM7091457, by = "GENE") %>% 
  full_join(GSM7091458, by = "GENE") %>% 
  tibble::column_to_rownames("GENE") %>% 
  as.matrix %>% 
  Matrix::Matrix(sparse = TRUE) %>% 
  CreateSeuratObject(project = "GSE227088")

scobj@meta.data$sample <- scobj@meta.data %>% rownames %>% str_remove_all("^[ACTG]*-")

table(scobj@meta.data$sample)
# 1-GSM7091453 2-GSM7091454 3-GSM7091455 4-GSM7091456 5-GSM7091457 6-GSM7091458
# 4576         5919         5531         6418         5184         5672
# QC----------------------------------------------------------------------------
# We applied the down sample analysis among samples sequenced according to the 
# mapped barcoded reads per cell of each sample and finally achieved the aggregated 
# matrix. Cells contained over 200 expressed genes and mitochondria UMI rate 
# below 10% passed the cell quality filtering and mitochondria genes were removed 
# in the expression table. 

scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^mt-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^Rp[sl]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^Hb[^(p)]")

pdf("QC.pdf")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nFeature = scobj@meta.data$nFeature_RNA,
  group = scobj@meta.data$orig.ident
) %>% 
  ggplot(aes(x = nFeature, color = group)) +
  geom_density() + 
  scale_x_continuous(breaks = c(200,300,400,500,1000,3000,4000,5000))
dev.off()

scobj <- subset(
  scobj, 
  subset = nFeature_RNA > 500 & 
    nFeature_RNA < 4000 & 
    # percent.mt < 10 & 
    percent.rp < 5)  

pdf("QC-2.pdf")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)
dev.off()

scobj <- scobj %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA(npcs = 50)
# Integration-------------------------------------------------------------------
# Perform integration (harmony)
scobj.harmony <- RunHarmony(
  scobj,
  group.by.vars = "sample",
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
DimPlot(scobj.harmony, group.by = "sample", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.5", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "sample", reduction = "tsne", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.5", reduction = "tsne", label = T) 
dev.off()
# Filter doublets---------------------------------------------------------------
if(TRUE){
  scobj.harmony.sce <- as.SingleCellExperiment(scobj.harmony)
  scobj.harmony.sce <- scDblFinder(scobj.harmony.sce, dbr.sd=1)
  
  # observed. We quantified the batch effect using Shannon entropy
  # as described in the CellMixS R package47. Mean entropy of the transcript data
  # was 0.91 and 0.92 for ADT data (Supplementary Fig. 4c,d). We observed that cells
  # from a single cluster (indicated by the black arrow) had lower entropy values and
  # contributed most to the batch effect.
  scobj.harmony.sce <- CellMixS::entropy(
    scobj.harmony.sce, 
    group = "sample", 
    k = 20, 
    res_name = "entropy_sample")
  scobj.harmony.sce <- CellMixS::entropy(
    scobj.harmony.sce, 
    group = "RNA_snn_res.0.5", 
    k = 20, 
    res_name = "entropy_RNA_snn_res.0.5")
  
  scobj.harmony <- as.Seurat(scobj.harmony.sce)
  scobj.h.sc <- subset(scobj.harmony, scDblFinder.class == "singlet")
  
  scobj.h.sc@meta.data %>%
    group_by(sample) %>%
    mutate(mean = mean(entropy_sample)) %>%
    distinct(mean)

  scobj.h.sc@meta.data %>%
    group_by(RNA_snn_res.0.5) %>%
    mutate(mean = mean(entropy_RNA_snn_res.0.5)) %>%
    distinct(mean)
}
# Annotation--------------------------------------------------------------------
# In our study, we found 28 clusters of sequenced cells in cardiac tissues and 
# identified all major cell types, including cardiomyocytes, fibroblasts, 
# endothelial cells, macrophages, smooth muscle cells, T cells, and B cells: 
# B cells (Marker: Cd79a, Ms4a1, and Cd19), 
# cardiomyocyte (Marker: Tnnc2, Ttn, and Pln), 
# endothelium (Marker: Pecam1, Cdh5, and Npr3), 
# fibroblast (Marker: Dcn, Fbln1, and Col1a2). 
# Myeloid cells, including macrophage (Marker: C1qb, Adgre1, and Ms4a7), 
# monocyte (Marker: Vcan, Ccr2, and Cd14), 
# neutrophil (Marker: Csf3r, S100a8, and Retnlg), 
# smooth muscle cells (Marker: Notch3, Pdgfrb, and Cspg4) and 
# T cells (Marker: Cd3e, Cd2, and Cd3d) (Figures 2C, D). 
# Gene expression analysis showed distinct patterns of high expression of the 14 
# transcripts in each of these cell types (Figure 2E). Our findings are consistent 
# with published studies on molecular profiles and cellular composition . 
# We have attempted to integrate SBP data to establish a new cell landscape that 
# might provide a new research strategy and evidence for the multi-target regulation of SBP.

scobj.h.sc <- RegroupIdents(scobj.h.sc, metadata = "RNA_snn_res.0.5")
scobj.h.sc -> scobj.h.sc.copy

CM <- c('Pln', 'Ttn', 'Tnnc2') # 1 9
ACM <- c('Tnnt2', 'Ttn', 'Myh6', 'Myh7') # 1 9
VCM <- c('Myh7', 'Myl2', 'Fhl2') # 1 9
EC <- c('Vwf', 'Cdh5', 'Thbd', 'Tek', 'Eng', 'Pecam1', 'Npr3') # 3 4 ?6 ?8
FB <- c('Fn1', 'Fbln1', 'Vim', 'Col1a1', 'Col1a2', 'Dcn') # 5
M <- c('Cd163', 'S100a8', 'Csf1r', 'C5ar1', 'Cd74') # 2
Myeloid	<- c('Cd163', 'Cd14', 'C1qa', 'C1qb', 'Adgre1', 'Ms4a7') # 2
Monocyte <- c('Vcan', 'Ccr2', 'Cd14') # 2
neutrophil <- c('Csf3r', 'S100a8', 'Retnlg') # 12
P <- c('Hcn1', 'Hcn4', 'Cacna1d') # pacemaker cells (P cells)
Peri <- c('Rgs5', 'Abcc9', 'Kcnj8') # 
SMC <- c('Myh11', 'Acta2','Cnn1', 'Tagln', 'Myocd', 'Cald1', 'Mylk', 'Notch3', 
         'Pdgfrb', 'Cspg4') #
Meso <- c('Msln', 'Wt1', 'Bnc1') # 
Neural <- c('Plp1', 'Nrxn1', 'Nrxn3') 
Adi <- c('Gpam', 'Fasn', 'Lep') # 4
Lymph <- c('Cd3e', 'Il7r', 'Cd40lg') # 7
T_cell <- c('Cd3e', 'Cd3d', 'Cd2', 'Gzmk') # 7
B_cell <- c('Cd79a', 'Ms4a1', 'Cd19') # 11
NK <- c('Gzmb', 'Rpf1') # 7
Mast <- c('Kit', 'Cpa3') 
RBC <- c('Hba1', 'Hba2', 'Hbb', 'Hbd', 'Hbe1', 'Hbg1', 'Hbg2', 'Hbm', 'Hbq1', 'Hbz') # 
CellCycle	<- c('Mki67', 'Ccna2', 'Ccnb2', 'Pcna', 'Stmn1')

VlnPlot(scobj.h.sc, features = CellCycle, group.by = "RNA_snn_res.0.5", layer = "counts", log = TRUE)

scobj.h.sc@meta.data$cell_type <- "Unknown"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.5 %in% c(1,9)] <- "CM"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.5 %in% c(2)] <- "Myeloid" 
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.5 %in% c(3)] <- "EC" 
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.5 %in% c(4, 6, 8)] <- "EC" 
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.5 %in% c(0, 5)] <- "FB" 
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.5 %in% c(7)] <- "Lymphoid" 
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.5 %in% c(11)] <- "B"  
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$RNA_snn_res.0.5 %in% c(12)]<- "Neutrophil" 

scobj.h.sc <- subset(scobj.h.sc, subset = cell_type != "Unknown")

p1 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "RNA_snn_res.0.5",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "sample",
              label = TRUE, pt.size = 0.5) 
p1 + p2 + p3

save(scobj.h.sc, file = "GSE227088-scobj.h.sc.annot.Rdata")




