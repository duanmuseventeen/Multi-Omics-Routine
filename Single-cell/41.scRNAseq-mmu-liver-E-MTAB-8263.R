# Acute liver failure is regulated by MYC- and microbiome-dependent programs---

# Data availability
# The small cytoplasmic RNA-seq data have been deposited with ArrayExpress
# under accession no. E-MTAB-8263, and 16S sequencing data with the European
# Nucleotide Archive under accession no. ERP116956. Source data are provided
# with this paper.

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
# E-MTAB-8263-------------------------------------------------------------------
scobj <- "scRNAseq" %>% 
  Read10X %>% 
  CreateSeuratObject(assay = "RNA", project = "E-MTAB-8263")

# scobj
# # An object of class Seurat 
# # 27998 features across 56527 samples within 1 assay 
# # Active assay: RNA (27998 features, 0 variable features)
# # 1 layer present: counts
# 
# table(str_extract(rownames(scobj@meta.data),"[0-9]*$"))
# # 1   10   11   12   13   14   15   16   17   18   19    2   20   21   22   23   24 
# # 3387 2137 1611  300 2520 2780 3058 2251  666  871 2048  358 1866 2015 1464 2360 2661 
# # 25   26   27   28   29    3   30    4    5    6    7    8    9 
# # 870 1446 1431 2213 1752 2517 2042 3093  509 2990 1077 1919 2315 
# add meta.data-----------------------------------------------------------------
meta <- data.table::fread("aggregation_sample_order.csv") 

meta.data <- scobj@meta.data
meta.data$barcode <- meta.data %>% rownames
meta.data <- meta.data %>% 
  mutate(order = str_extract(barcode, "[0-9]*$")) %>% 
  left_join(meta %>% mutate(order = as.character(order)), by = "order")
rownames(meta.data) <- meta.data$barcode

scobj@meta.data <- meta.data

# table(scobj@meta.data$library_id)
# # ABX_APAP_1        ABX_APAP_2              GF_1              GF_2 
# # 1611               300              3387               358 
# # GF_3         GF_APAP_1         GF_APAP_2          GF_TAA_1 
# # 2517              2315              2137              2048 
# # GF_TAA_2      MyD88_Trif_1      MyD88_Trif_2 MyD88_Trif_APAP_1 
# # 1866              1431              2213              1752 
# # MyD88_Trif_APAP_2             SPF_1             SPF_2             SPF_3 
# # 2042              3093               509              2990 
# # SPF_APAP_1        SPF_APAP_2        SPF_APAP_3        SPF_APAP_4 
# # 2520              2780              3058              2251 
# # SPF_APAP_KJPyr9_1 SPF_APAP_KJPyr9_2      SPF_KJPyr9_1      SPF_KJPyr9_2 
# # 666               871              1077              1919 
# # SPF_TAA_1         SPF_TAA_2         SPF_TAA_3         SPF_TAA_4 
# # 2015              1464              2360              2661 
# # SPF_TAA_KJPyr9_1  SPF_TAA_KJPyr9_2 
# # 870              1446 
# QC----------------------------------------------------------------------------
scobj <- subset(
  scobj, 
  subset = library_id %in% c("SPF_1","SPF_2","SPF_3","SPF_APAP_1","SPF_APAP_2",
                             "SPF_APAP_3","SPF_APAP_4","SPF_TAA_1","SPF_TAA_2",
                             "SPF_TAA_3","SPF_TAA_4"))

table(scobj@meta.data$library_id)
# SPF_1      SPF_2      SPF_3 SPF_APAP_1 SPF_APAP_2 SPF_APAP_3 SPF_APAP_4  SPF_TAA_1 
# 3093        509       2990       2520       2780       3058       2251       2015 
# SPF_TAA_2  SPF_TAA_3  SPF_TAA_4 
# 1464       2360       2661 

scobj <- subset(
  scobj, 
  subset = library_id %in% c("SPF_1","SPF_3","SPF_APAP_1","SPF_APAP_2",
                             "SPF_APAP_3","SPF_APAP_4","SPF_TAA_1","SPF_TAA_2",
                             "SPF_TAA_3","SPF_TAA_4"))

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$library_id)

scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^mt-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^Rp[sl]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^Hb[^(p)]")

# VlnPlot(scobj,
#         group.by = "library_id",
#         features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
# FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# data.frame(
#   nFeature = scobj@meta.data$nFeature_RNA,
#   group = scobj@meta.data$orig.ident
# ) %>% 
#   ggplot(aes(x = nFeature, color = group)) +
#   geom_density() + 
#   scale_x_continuous(breaks = c(200,300,400,500,1000,3000,4000,5000))

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 200 & 
    nFeature_RNA < 2000 & 
    nCount_RNA < 5000 &
    percent.mt < 10 &
    percent.hb < 1)  

# scobj.harmony@assays$RNA$counts.SPF_1
# ...
# scobj.harmony@assays$RNA$counts.SPF_TAA_4
# 
# scobj.harmony@assays$SCT$counts
# scobj.harmony@assays$SCT$data
# scobj.harmony@assays$SCT$scale.data

# VlnPlot(scobj,
#         group.by = "library_id",
#         features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
        
scobj.harmony <- scobj %>% 
  SCTransform %>% # vst.flavor = 'v2', verbose = FALSE
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = "library_id",
    reduction.use = "pca",
    reduction.save = "harmony") 

scobj.harmony <- scobj.harmony %>% 
  JoinLayers(assay = "RNA")

scobj.harmony@assays$SCT$counts[c("Sox17","Mrpl15"),1:20]
# 2 x 20 sparse Matrix of class "dgCMatrix"
# [[ suppressing 20 column names ‘AAACCTGAGCAGCCTC-4’, ‘AAACCTGAGGGAACGG-4’, ‘AAACCTGGTAACGCGA-4’ ... ]]
# 
# Sox17  . . . . . . . . . . . . . . . . . . . .
# Mrpl15 . 2 . . . . 1 1 . . . . . . . . . 1 . .
scobj.harmony@assays$RNA$counts[c("Sox17","Mrpl15"),1:20]
# 2 x 20 sparse Matrix of class "dgCMatrix"
# [[ suppressing 20 column names ‘AAACCTGAGCAGCCTC-4’, ‘AAACCTGAGGGAACGG-4’, ‘AAACCTGGTAACGCGA-4’ ... ]]
# 
# Sox17  . . . . . . . . . . . . . . . . . . . .
# Mrpl15 . 1 . . . . 1 1 . . . . . . . . . 1 . .

scobj.harmony <- scobj.harmony %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1, 1, 0.1)) %>% 
  RunUMAP(reduction = "harmony", min_dist = 0.3, dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30)

clustree(scobj.harmony, prefix = "SCT_snn_res.")

scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "SCT_snn_res.0.4")
DimPlot(scobj.harmony, group.by = "library_id", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "SCT_snn_res.0.4", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "SCT_snn_res.0.4", reduction = "tsne", label = T)
# Cell cycle--------------------------------------------------------------------
# annot <- AnnotationDbi::select(
#   org.Mm.eg.db, 
#   keys = keys(org.Mm.eg.db,keytype = "ENSEMBL"),
#   columns = "SYMBOL", 
#   keytype = "ENSEMBL")
# 
# dat <- data.frame(
#   SYMBOL = rownames(scobj.harmony),
#   stringsAsFactors = FALSE
# ) %>% 
#   left_join(annot, by = "SYMBOL")
# 
# sum(is.na(dat$ENSEMBL))
# # [1] 1465
# 
# dat[is.na(dat$ENSEMBL),] %>% head

annot <- data.table::fread("scRNAseq/genes.tsv", col.names = c("ENSEMBL","SYMBOL","TYPE"))

dat <- data.frame(
  SYMBOL = rownames(scobj.harmony),
  stringsAsFactors = FALSE
) %>%
  left_join(annot %>% dplyr::select(-TYPE), by = "SYMBOL")

sum(is.na(dat$ENSEMBL))
# [1] 7

library(scran)
# hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
str(mm.pairs)
head(mm.pairs$G1)

all(rownames(scobj.harmony) == dat$SYMBOL)

# sce <- as.SingleCellExperiment(scobj.harmony, assay = "RNA")
# 警告信息:
# 1: Layer ‘data’ is empty 
# 2: Layer ‘scale.data’ is empty 

# assignments <- cyclone(sce, pairs = mm.pairs)

# head(assignments$scores)

# table(assignments$phases)

# doublets----------------------------------------------------------------------
if(TRUE){
  # https://www.jianshu.com/p/21b2f08652ed
  # https://blog.csdn.net/hx2024/article/details/140810266
  
  # scDblFinder
  sce.harmony <- as.SingleCellExperiment(scobj.harmony, assay = "SCT")
  sce.harmony <- scDblFinder(
    sce.harmony,
    samples = sce.harmony$library_id,
    clusters = sce.harmony$SCT_snn_res.0.4,
    dbr.sd = 1) 
  
  scobj.harmony[['scDblFinder.class']] <- sce.harmony$scDblFinder.class
  
  # table(scobj.harmony@meta.data$scDblFinder.class) # use "RNA" layers
  # singlet doublet 
  # 21194     585 
  
  table(scobj.harmony@meta.data$scDblFinder.class) # use "SCT" layers
  # singlet doublet 
  # 21648     131 
}
# Re-run------------------------------------------------------------------------
scobj[['scDblFinder.class']] <- sce.harmony$scDblFinder.class

scobj.harmony <- scobj %>% 
  subset(subset = scDblFinder.class == "singlet") %>%
  SCTransform %>% # vst.flavor = 'v2', verbose = FALSE
  RunPCA(npcs = 50) %>% 
  RunHarmony(
    group.by.vars = "library_id",
    reduction.use = "pca",
    reduction.save = "harmony") %>% 
  JoinLayers(assay = "RNA") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1, 1, 0.1)) %>% 
  RunUMAP(reduction = "harmony", min_dist = 0.3, dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30)

clustree(scobj.harmony, prefix = "SCT_snn_res.")

scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "SCT_snn_res.0.1")
DimPlot(scobj.harmony, group.by = "library_id", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "SCT_snn_res.0.1", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "SCT_snn_res.0.1", reduction = "tsne", label = T)
# Annotation--------------------------------------------------------------------
# MP, mononuclear phagocytes
MP <-	c('Cd68',	'Itgam',	'Itgax',	'HLA+AC0-DRA',	'Csf1r',	'Cd14')	# 
Macrophage <- c("Cd68", "Cd163", "Lyz2") # 
Myeloid	<- c('Cd163', 'Cd14', 'C1qa', 'C1qb', 'Adgre1', 'Ms4a7') # 
Monocyte <- c('Vcan', 'Ccr2', 'Cd14') # 2
neutrophil <- c('Csf3r', 'S100a8', 'Retnlg') # 
B_cell <- c('Cd79a', 'Ms4a1', 'Cd19') # 
Mast <- c('Kit', 'Cpa3') 
hept <- c("Cyp3a7","Cyp2a7","Cyp2a6","Scd","Hmgcs1","Acss2",
          "Slbp","Rnd3","Pck1","G6pc","Ghr","Aldh6a1","Rpp25l","Hsd11b1",
          "Hamp","Ghr","Hpr","Akr1c1","Masp2") # 
EC <- c('Vwf', 'Cdh5', 'Thbd', 'Tek', 'Eng', 'Pecam1', 'Npr3') #
chol <- c("Krt7","Krt19","Sox9","Epcam") # 
DC <- c('Klrd1', 'Flt3', 'H2-Ab1') # 9
stel <- c("Acta2","Col1a1","Rbp1","Rgs5","Lrat","Pdgfrb","Pdgfra") # 0 8
T_cell <- c('Cd3e', 'Cd3d', 'Cd2', 'Gzmk') # 
fibroblasts <- c("Fn1",'Dcn', 'Sfrp4', 'Lum', 'Col1a1','Col1a1', 'Pclaf') # 
# CellCycle	<- c('Mki67', 'Ccna2', 'Ccnb2', 'Pcna', 'Stmn1')

VlnPlot(scobj.harmony, features = CellCycle, group.by = "SCT_snn_res.0.1", layer = "counts", log = TRUE)

scobj.harmony@meta.data$cell_type <- "Unknown"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(0)] <- "Stellate"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(1)] <- "Endothelial"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(2)] <- "Monocyte"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(3)] <- "Neutrophil"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(4)] <- "T"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(5)] <- "Macrophage"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(6)] <- "B"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(7)] <- "Hepatocyte"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(8)] <- "Stellate"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(9)] <- "Dendritic cell"
scobj.harmony@meta.data$cell_type[scobj.harmony@meta.data$SCT_snn_res.0.1 %in% c(10)] <- "Cholangiocyte"

p1 <- DimPlot(scobj.harmony, reduction = "umap",
              group.by = "SCT_snn_res.0.1",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj.harmony, reduction = "umap",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj.harmony, reduction = "umap",
              group.by = "library_id",
              label = TRUE, pt.size = 0.5) 
p4 <- FeaturePlot(scobj.harmony, reduction = "umap", features = "B4galt1")
p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2)

table(scobj.harmony@meta.data$cell_type)
# B  Cholangiocyte Dendritic cell    Endothelial     Hepatocyte 
# 993             62            118           4435            283 
# Macrophage       Monocyte     Neutrophil       Stellate              T 
# 1488           3535           3031           6059           1644

RidgePlot(object = scobj.harmony, group.by = "cell_type", features = 'B4galt1')
VlnPlot(object = scobj.harmony, group.by = "cell_type", features = 'B4galt1')
DotPlot(scobj.harmony, group.by = "cell_type", features = 'B4galt1')

save(scobj.harmony, file = "E-MTAB-8263-scobj.harmony.annot.Rdata")


