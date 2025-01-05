# Single-Cell Transcriptome Analysis Reveals Dynamic Cell Populations and Differential Gene Expression Patterns in Control and Aneurysmal Human Aortic Tissue---- 

# GSE155468
# BACKGROUND: Ascending thoracic aortic aneurysm (ATAA) is caused by the progressive weakening and dilatation of the aortic wall and can lead to aortic dissection, rupture, and other life-threatening complications. To improve our understanding of ATAA pathogenesis, we sought to comprehensively characterize the cellular composition of the ascending aortic wall and to identify molecular alterations in each cell population of human ATAA tissues.
# METHODS: We performed single-cell RNA sequencing (sc-RNAseq) analysis of ascending aortic tissues from 11 study participants, including 8 patients with ATAA (4 women and 4 men) and 3 controls (2 women and 1 man). Cells extracted from aortic tissue were analyzed and categorized by using sc-RNAseq data to perform cluster identification. ATAA-related changes were then examined by comparing the proportions of each cell type and the gene expression profiles between ATAA and control tissues. We also examined which genes may be critical for ATAA by performing the integrative analysis of our sc-RNAseq data with publicly available data from genome-wide association studies (GWAS).
# RESULTS: We identified 11 major cell types in human ascending aortic tissue; the high-resolution reclustering of these cells further divided them into 40 subtypes. Multiple subtypes were observed for smooth muscle cells, macrophages, and T lymphocytes, suggesting that these cells have multiple functional populations in the aortic wall. Generally, ATAA tissues had fewer nonimmune cells and more immune cells, especially T lymphocytes, than did control tissues. Differential gene expression data suggested the presence of extensive mitochondrial dysfunction in ATAA tissues. In addition, integrative analysis of our sc-RNAseq data with public GWAS data and promoter capture Hi-C data suggested that ERG (ETS [erythroblast transformation-specific] related gene) exerts an important role in maintaining normal aortic wall function.
# CONCLUSIONS: Our study provides a comprehensive evaluation of the cellular composition of the ascending aortic wall and reveals how the gene expression landscape is altered in human ATAA tissue. The information from this study makes important contributions to our understanding of ATAA formation and progression.
# mRNA profiles of 3 control and 8 ATAA human aortic tissue at single cell level

Variable	ATAA1	ATAA2	ATAA3	ATAA4	ATAA5	ATAA6	ATAA7	ATAA8	Control4	Control6	Control9
Sex	F	F	M	M	F	F	M	M	F	M	F
Ethnicity	Non-Hispanic	Non-Hispanic	Non-Hispanic	Non-Hispanic	Non-Hispanic	Non-Hispanic	Non-Hispanic	Non-Hispanic	Non-Hispanic	Non-Hispanic	Hispanic
Race	White	White	White	White	White	White	White	White	White	Black	Latino
Age (y)	75	78	59	62	75	67	69	56	63	61	62
Diagnosis/ Comments	ATAA	ATAA	ATAA with root aneurysm	ATAA	ATAA with root aneurysm	ATAA with arch and DTAA	ATAA with root aneurysm	ATAA with root aneurysm	Heart transplant recipient	Heart transplant recipient	Lung transplant donor
Aortic diameter (cm)	5.2	4.9	5	5.2	5.8	4.9	5.2	5.2	NA	NA	2.2
Smoking status	Past (quit before 1990)	Never	Never	Never	Past (quit 1999)	Never	Never	Never	Never	Past	Current
Diabetes	No	No	Yes	Yes	No	No	No	No	No	No	Yes
Hypertension	Yes	Yes	Yes	Yes	Yes	Yes	Yes	Yes	No	Yes	Yes
COPD	No	Yes	No	No	Yes	No	No	No	No	No	No
Aortic valve regurgitation	No	No	Yes	Yes	Yes	Yes	No	Yes	No	No	No
BAV	No	Yes	Yes	No	No	No	NA	No	NA	NA	No
Re-operation	No	No	Yes *	No	No	No	Yes **	No	No	No	No

# Reference---------------------------------------------------------------------
# https://cloud.tencent.com/developer/article/2317475
# https://www.jianshu.com/p/784a04c49873
# https://blog.csdn.net/weixin_40695088/article/details/136358659
# https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_GLYCOSYLATION.html

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
hub <- dir()
names(hub) <- hub

scobj.list <- lapply(hub, function(x){
  dat <- data.table::fread(x) %>% 
    tibble::column_to_rownames("V1")
  colnames(dat) <- colnames(dat) %>% paste0("_",x) %>% str_remove(".txt.gz$")
  dat %>% as.matrix %>% 
    Matrix::Matrix(sparse = TRUE) %>%
    CreateSeuratObject(counts = ., project = x %>% str_remove(".txt.gz$")) %>% 
    return
})

scobj.list

scobj <- merge(x = scobj.list[[1]], y = scobj.list[-1])

scobj

table(scobj@meta.data$orig.ident)
# GSM4704931_Con4 GSM4704932_Con6 GSM4704933_Con9 GSM4704934_TAA1 GSM4704935_TAA2 
# 3377            1193            3907            6686            4598 
# GSM4704936_TAA3 GSM4704937_TAA4 GSM4704938_TAA5 GSM4704939_TAA6 GSM4704940_TAA7 
# 3622            6384            4802            3526            6997 
# GSM4704941_TAA8 
# 3036 
# QC----------------------------------------------------------------------------
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^HB[^(P)]")

# VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
# FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# data.frame(
#   nFeature = scobj@meta.data$nFeature_RNA,
#   group = scobj@meta.data$orig.ident
# ) %>%
#   ggplot(aes(x = nFeature, color = group)) +
#   geom_density() +
#   scale_x_continuous(breaks = c(200,300,400,500,750,1000,3000,4000,5000))

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 750 &
    nFeature_RNA < 5000 & 
    nCount_RNA < 20000 &
    percent.mt < 10 &
    percent.rp < 40 &
    percent.hb < 1)  

VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)

table(scobj@meta.data$orig.ident)
# GSM4704931_Con4 GSM4704932_Con6 GSM4704933_Con9 GSM4704934_TAA1 GSM4704935_TAA2 
# 2931            1006            3666            4988            3191 
# GSM4704936_TAA3 GSM4704937_TAA4 GSM4704938_TAA5 GSM4704939_TAA6 GSM4704940_TAA7 
# 3166            5470            4370            2754            6522 
# GSM4704941_TAA8 
# 2669 

scobj.ds <- subset(scobj, subset = orig.ident != "GSM4704932_Con6", downsample = 2600)

table(scobj.ds@meta.data$orig.ident)
# GSM4704931_Con4 GSM4704933_Con9 GSM4704934_TAA1 GSM4704935_TAA2 GSM4704936_TAA3 
# 2600            2600            2600            2600            2600 
# GSM4704937_TAA4 GSM4704938_TAA5 GSM4704939_TAA6 GSM4704940_TAA7 GSM4704941_TAA8 
# 2600            2600            2600            2600            2600 

scobj.ds <- scobj.ds %>% 
  SCTransform %>% # vst.flavor = 'v2', verbose = FALSE
  RunPCA(npcs = 50)
# Integration-------------------------------------------------------------------
scobj.ds.harmony <- RunHarmony(
  scobj.ds,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  reduction.save = "harmony")
# Cluster-----------------------------------------------------------------------
scobj.ds.harmony <- FindNeighbors(scobj.ds.harmony, reduction = "harmony", dims = 1:30)
scobj.ds.harmony <- FindClusters(scobj.ds.harmony, resolution = seq(0.1, 1, 0.1))

clustree(scobj.ds.harmony, prefix = "SCT_snn_res.")
# Visualization-----------------------------------------------------------------
#UMAP
scobj.ds.harmony <- RunUMAP(scobj.ds.harmony, reduction = "harmony", min_dist = 0.3, dims = 1:30)
#T-SNE
scobj.ds.harmony <- RunTSNE(scobj.ds.harmony, reduction = "harmony", dims = 1:30)

scobj.ds.harmony@meta.data <- scobj.ds.harmony@meta.data %>% 
  mutate(disease = str_extract(orig.ident, "[a-zA-Z0-9]+$") %>% str_remove("[0-9]"))

DimPlot(scobj.ds.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.ds.harmony, group.by = "disease", reduction = "umap", label = T)
DimPlot(scobj.ds.harmony, group.by = "SCT_snn_res.0.3", reduction = "umap", label = T)
DimPlot(scobj.ds.harmony, group.by = "SCT_snn_res.0.3", reduction = "tsne", label = T)
# Cell cycle--------------------------------------------------------------------
exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scobj.ds.harmony <- CellCycleScoring(scobj.ds.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(scobj.ds.harmony, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

scobj.ds.harmony <- RunPCA(scobj.ds.harmony, features = c(s.genes, g2m.genes))
DimPlot(scobj.ds.harmony, reduction = "pca")
DimPlot(scobj.ds.harmony, reduction = "umap")
# Filter doublets---------------------------------------------------------------
if(TRUE){
  # scDblFinder
  scobj.ds.harmony.sce <- as.SingleCellExperiment(scobj.ds.harmony, assay = "SCT")
  scobj.ds.harmony.sce <- scDblFinder(
    scobj.ds.harmony.sce, 
    samples = scobj.ds.harmony.sce$orig.ident,
    dbr.sd=1)
  
  scobj.ds.harmony <- as.Seurat(scobj.ds.harmony.sce)
  scobj.h.sc <- subset(scobj.ds.harmony, scDblFinder.class == "singlet")
}
# AUCELL------------------------------------------------------------------------
require(GSEABase)

#### INFLAMATION--------------------------------------------------------------
AIR <- getGmt("https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=GOBP_ACTIVATION_OF_IMMUNE_RESPONSE&fileType=gmt")

# AIR@.Data[[1]]@geneIds %>% length
# data.table::fread("GOBP_ACTIVATION_OF_IMMUNE_RESPONSE.v2024.1.Hs.gmt") %>% 
#   colnames %>% 
#   length

# select counts as input
cell_rankings <- AUCell_buildRankings(
  scobj.h.sc@assays$SCT$counts,
  plotStats = TRUE,
  verbose = TRUE)
# Quantiles for the number of genes detected by cell: 
#   (Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
# min   1%   5%  10%  50% 100% 
# 208  780  876  959 1618 4540 

# 1
cell_AUC <- AUCell_calcAUC(geneSets = AIR, rankings = cell_rankings, aucMaxRank = 5)
# 2
genelist <- data.table::fread("GOBP_ACTIVATION_OF_IMMUNE_RESPONSE.v2024.1.Hs.gmt") %>% 
  colnames
cell_AUC <- AUCell_calcAUC(
  geneSets = genelist[-c(1:2)], 
  rankings = cell_rankings, 
  aucMaxRank = ceiling(0.05 * nrow(cell_rankings)))

# geneSets: List of gene-sets (or signatures) to test in the cells. The gene-sets should be provided as GeneSet, GeneSetCollection or character list (see examples).
# aucMaxRank: In a simplified way, the AUC value represents the fraction of genes, within the top X genes in the ranking, that are included in the signature. The parameter 'aucMaxRank' allows to modify the number of genes (maximum ranking) that is used to perform this computation. By default, it is set to 5% of the total number of genes in the rankings. Common values may range from 1 to 20%.

# Genes in the gene sets NOT available in the dataset: 
#   GOBP_ACTIVATION_OF_IMMUNE_RESPONSE: 	63 (11% of 583)

AUCell_plotHist(cell_AUC["GOBP_ACTIVATION_OF_IMMUNE_RESPONSE",], nBreaks = 40)

cell_assignment <- AUCell_exploreThresholds(cell_AUC, plotHist = TRUE, assign=TRUE)

# newAssignments <- AUCell_assignCells(cell_AUC, getThresholdSelected(cell_assignment))
# getAssignments(newAssignments)

# getThresholdSelected is equal to cell_assignment$GOBP_ACTIVATION_OF_IMMUNE_RESPONSE$aucThr$selected

scobj.h.sc@meta.data$id <- rownames(scobj.h.sc@meta.data)
scobj.h.sc@meta.data <- scobj.h.sc@meta.data %>% 
  mutate(GOBP_ACTIVATION_OF_IMMUNE_RESPONSE = case_when(
    id %in% cell_assignment$GOBP_ACTIVATION_OF_IMMUNE_RESPONSE$assignment ~ TRUE,
    TRUE ~ FALSE))

DimPlot(scobj.h.sc, group.by = "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE", label = TRUE)
#### GLYCOSYLATION--------------------------------------------------------------
glycosylation <- getGmt("https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=GOBP_GLYCOSYLATION&fileType=gmt")

# select counts as input
cell_rankings <- AUCell_buildRankings(scobj.h.sc@assays$SCT$counts)
# Quantiles for the number of genes detected by cell: 
#   (Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
# min   1%   5%  10%  50% 100% 
# 208  780  876  959 1618 4540 

cell_AUC <- AUCell_calcAUC(geneSets = glycosylation, rankings = cell_rankings)
# Genes in the gene sets NOT available in the dataset: 
#   GOBP_GLYCOSYLATION: 	28 (11% of 244)

cell_assignment <- AUCell_exploreThresholds(cell_AUC, plotHist = TRUE, assign=TRUE)

cell_assignment$GOBP_GLYCOSYLATION$aucThr$selected

scobj.h.sc@meta.data$id <- rownames(scobj.h.sc@meta.data)
scobj.h.sc@meta.data <- scobj.h.sc@meta.data %>% 
  mutate(GO_GLYCOSYLATION = case_when(
    id %in% cell_assignment$GOBP_GLYCOSYLATION$assignment ~ TRUE,
    TRUE ~ FALSE))

DimPlot(scobj.h.sc, group.by = "GO_GLYCOSYLATION", label = TRUE)
# AddModuleScore----------------------------------------------------------------
require(msigdbr)

AIR <- getGmt("https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=GOBP_ACTIVATION_OF_IMMUNE_RESPONSE&fileType=gmt")
AIR@.Data[[1]]@geneIds

# if anyone gene in gmt are not involved in scobj, the function cannot be run
AIR.filtered <- AIR@.Data[[1]]@geneIds[AIR@.Data[[1]]@geneIds %in% rownames(scobj.h.sc)]

scobj.h.sc <- AddModuleScore(
  scobj.h.sc,
  features = list(AIR.filtered), # the format of input geneset is vital
  ctrl = 100,
  name = "SEURAT_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE")

scobj.h.sc@meta.data$SEURAT_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE1

FeaturePlot(scobj.h.sc, features = "SEURAT_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE1")
# Annotation--------------------------------------------------------------------
scobj.h.sc <- RegroupIdents(scobj.h.sc, metadata = "SCT_snn_res.0.3")

# The most common cell type observed in aortas were VSMCs, which comprised 
# approximately ~65% of nuclei in our samples (clusters 1, 2, 3), consistent 
# with observations from aortic histology. (Figure 2c, Supplemental Table 3) 
# Clusters 1, 2, 3, and 4 all expressed canonical markers of VSMCs including 
# MYH11 and ACTA2, and included VSMC subtypes as well as pericytes.

# Fibroblast cluster 6 was the next most common cell type after VSMCs and expressed typical but less specific gene markers such as LAMA2 and TSHZ2.

# As expected for vascular tissue, endothelial cells expressing PECAM1 (CD31) were also numerous and were observed in three separate clusters 8, 9, 10.

SMC <- c('MYH11', 'ACTA2','CNN1', 'SMTN', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') # 4 ?5 7 9 ?15
# Peri <- c('RGS5', 'ABCC9', 'KCNJ8') # 4 5 7 9
FB <- c('LAMA2', 'TSHZ2', "FN1", "COL1A1") # 7 ?9
EC <- c('PECAM1',"VWF") # 11
RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ') # 
Meso <- c('MSLN', 'WT1', 'BNC1') # 
Adi <- c('GPAM', 'FASN', 'LEP') # 
Neural <- c('PLP1', 'NRXN1', 'NRXN3') # 
Lymph <- c('CD3E', 'IL7R', 'CD40LG') # 0 2 ?6 10 12 13 15 16
Myeloid <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74', 'CD14', 'C1QA') # 1 3
Mast <- c("GATA2","HPGDS","HPGD","KIT","IL1RL1","LAPTM4A","MLPH","VWA5A","RGS13","LTC4S") # 14

if(TRUE){
  T_cell <- c("PTPRC","CD3D","IL7R","GZMK","IL32","CXCR4","CD52","TRAC","CD2","TRBC2")
  MonoMaphDC <- c("LYZ","CXCL8","IL1B","HLA-DRA","C1QA","C1QB","CXCL3","HLA-DRA1",
                  "CXCL2","C1QC","CD74","SELENOP")
  SMC <- c("ACTA2","MYL9","TAGLN","TPM2","CALD1","MYH11","DSTN","PLN","MFGE8",
           "MYH10","MAP1B","IGFBP2")
  NK <- c("XCL2","XCL1","KLRD1","GNLY","KLRC1","CD7","CTSW","PRF1","CMC1")
  Fibroblast <- c("DCN","LUM","FN1","MGP","CFH","COL1A2","AEBP1","IGFBP6","GSN","C1R")
  EC <- c("VWF","IF127","PECAM1","AQP1","SRRY1","SPARCL1","GNG11","IGFBP4","CAVIN2","ADGRL4")
  # mesenchymal stem cells
  MSC <- c("C11orf96","IGFBP5","CRISPLD2","MT2A","ADIRF","MT1M","PLAC9","PLAC9","NDUFA4L2")
  # plasma <- c("IGHA1","MZB1","SSR4","PRDX4","FKBP11","SEC11C","QPCT","LMF1")
  # B <- c("CD19","MS4A1","LY9","BANK1","ARHGAP24","CD79B","LTB","LINC00926","STAG3")
  Mast <- c("GATA2","HPGDS","HPGD","KIT","IL1RL1","LAPTM4A","MLPH","VWA5A","RGS13","LTC4S") # 14
}

VlnPlot(scobj.h.sc, features = NK, group.by = "SCT_snn_res.0.3", layer = "counts", log = TRUE)

scobj.h.sc@meta.data$cell_type <- "Unknown"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.3 %in% c(0,2,10,12,13,15,16)] <- "T cell"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.3 %in% c(1,3)] <- "Myeloid"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.3 %in% c(4,5,9)] <- "SMC"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.3 %in% c(6, 8)] <- "NK"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.3 %in% c(7)] <- "FB"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.3 %in% c(11)] <- "EC"
scobj.h.sc@meta.data$cell_type[scobj.h.sc@meta.data$SCT_snn_res.0.3 %in% c(14)] <- "Mast"

p1 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "SCT_snn_res.0.3",
              label = TRUE, pt.size = 0.5) 
p2 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "cell_type",
              label = TRUE, pt.size = 0.5)
p3 <- DimPlot(scobj.h.sc, reduction = "UMAP",
              group.by = "disease",
              label = TRUE, pt.size = 0.5) 
p1 + p2 + p3
# ggsave(..., width = 18, height = 6)

FeaturePlot(object = scobj.h.sc, slot = "counts", features = 'Targetgene')
RidgePlot(object = scobj.h.sc, features = 'Targetgene')
DotPlot(scobj.h.sc, group.by = "disease", features = "Targetgene")
DotPlot(scobj.h.sc, group.by = "cell_type", features = "Targetgene")

save(scobj.h.sc, file = "GSE155468-scobj.annot.Rdata")


