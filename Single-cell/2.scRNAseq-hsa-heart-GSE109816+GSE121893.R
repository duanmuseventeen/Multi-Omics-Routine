# Single-cell reconstruction of the adult human heart during heart failure and recovery reveals the cellular landscape underlying cardiac function
# https://www.nature.com/articles/s41556-019-0446-7#Sec9
# GSE109816 and GSE121893

# Samples from human participants
# This study was approved by the Ethics Committee of Fuwai Hospital, Chinese Academy of Medical Sciences and Peking Union Medical University, as well as by the Ministry of Science and Technology. Written informed consent for tissue donation, which clearly stated the purpose of our study, was obtained from all of the patients.
# 
# LV and/or LA tissue samples were collected from 14 previously healthy organ donors (12 males and 2 females). Right ventricular tissue samples were obtained from two additional previously healthy organ donors (male). All of the donors were deceased due to acute trauma or anoxia, and all of the hearts had a LV ejection fraction of greater than 50%, a cold ischaemia time of less than 6 h and were deemed to be not suitable for transplantation for technical or non-cardiac reasons. Atrial tissue samples (n = 9) were obtained from the human LA, whereas all of the ventricular tissue samples (n = 7) were dissected from the LV anterior wall close to the cardiac apex. The range of donor ages was 21–52 yr, with a median age of 45.5 yr. Information regarding the cells isolated from these tissue samples is provided in in Supplementary Table 1.
# 
# HF samples were collected from patients (6 males) undergoing heart transplantation. Atrial and ventricular tissue samples were dissected from the same regions as in normal hearts. Six cases of HF were included in our study, with two progressing from cHF and four progressing from dHF.
# 
# For heart samples from patients with HF who were treated with LVAD, two male patients with HF received continuous-flow LVAD implantation either as therapy or as a bridge to heart transplantation. The samples that were obtained before therapy were collected during the implantation procedure, which involved introducing an opening in the LV. Post-treatment samples were acquired during the removal of the device. The LV ejection fraction was less than 30% in both of the patients. One patient showed cardiac function recovery 166 d after LVAD implanting (LVAD1), whereas the other patient underwent cardiac transplantation 155 d after LVAD implantation owing to a lack of any signs of functional improvement (LVAD2). Ventricular tissue samples were collected during surgical operation of LVAD treatment.

# Note
# This datased is low-qualitied because of high percent of mito, low depth and 
# relatively higher RBC. I don't suggest use this datases to investigate.

# Reference---------------------------------------------------------------------
# https://lishensuo.github.io/posts/bioinfo/028%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--%E5%9F%BA%E4%BA%8E%E6%96%87%E7%8C%AE%E7%9A%84%E7%BB%86%E8%83%9E%E7%B1%BB%E5%9E%8B%E6%B3%A8%E9%87%8Amarker/
# http://117.50.127.228/CellMarker/
# http://www.360doc.com/content/12/0121/07/76149697_996050517.shtml
# https://cloud.tencent.com/developer/article/2401058

rm(list = ls());gc()

require(dplyr)
require(stringr)
require(Seurat)
require(SeuratObject)
require(Matrix)

# Load data---------------------------------------------------------------------
dat1 <- data.table::fread("GSE121893/GSE121893_human_heart_sc_umi.csv.gz")
dat2 <- data.table::fread("GSE121893/GSE109816_normal_heart_umi_matrix.csv.gz")
dat  <- dat1 %>% 
  left_join(dat2, by = "V1")

dim(dat1)
# [1] 25742  4934
dim(dat2)
# [1] 54750  9995
dim(dat)
# [1] 25742 14928

raw.data <- dat %>%
  tibble::column_to_rownames("V1") %>%
  as.data.frame

annot1 <- data.table::fread("GSE121893/GSE121893_human_heart_sc_info.txt.gz")
annot2 <- data.table::fread("GSE121893/GSE109816_normal_heart_cell_info.txt.gz")
annot <- bind_rows(annot1[,-35], annot2) %>% 
  mutate(orig.ident = ID, 
         disease = str_extract(Type, "^[HFN]{1,2}"),
         tissue = str_extract(Type, "[LAV]{1,2}"),
         CM = str_extract(Type, "[NCM]{2,3}$"))

all(annot$ID == colnames(raw.data))
# [1] TRUE
# Create Seurat Object---------------------------------------------------------
raw.data <- as.matrix(raw.data)
raw.data <- Matrix(raw.data, sparse = TRUE)
sce.all <- CreateSeuratObject(counts = raw.data, data = raw.data)
sce.all <- AddMetaData(object = sce.all, metadata = annot)
# QC-----------------------------------------------------------------------------
# sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all[["percent.rp"]] <- PercentageFeatureSet(sce.all, pattern = "^RP[SL]")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "mito.perc", "percent.rp"), ncol = 4)
FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

sce.all <- subset(sce.all, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.rp < 10)  
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "mito.perc", "percent.rp"), ncol = 4)

sce.all <- sce.all %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData %>% 
  RunPCA

#Examine and visualize PCA results a few different ways
print(sce.all[["pca"]], dims = 1:10, nfeatures = 5)
VizDimLoadings(sce.all, dims = 1:2, reduction = "pca")
DimPlot(sce.all, reduction = "pca")
DimHeatmap(sce.all, dims = 1:6, cells = 500, balanced = TRUE, )

# sce.all <- JackStraw(sce.all, dims = 50, num.replicate = 100)
# sce.all <- ScoreJackStraw(sce.all, dims = 1:50)
# JackStrawPlot(sce.all, dims = 1:50)
# ElbowPlot(sce.all, ndims = 50)

# Integrate Data----------------------------------------------------------------
sce.all.harmony <- RunHarmony(sce.all, 
                              group.by.vars = "Individual", 
                              reduction.use = "pca",
                              reduction.save = "harmony")

# Cluster-----------------------------------------------------------------------
sce.all.harmony <- FindNeighbors(sce.all.harmony, reduction = "harmony", dims = 1:25)
sce.all.harmony <- FindClusters(sce.all.harmony, resolution = seq(0.1, 1, 0.1)) 

clustree(sce.all.harmony, prefix = "RNA_snn_res.")
# Visualization-----------------------------------------------------------------
#UMAP
sce.all.harmony <- RunUMAP(sce.all.harmony, reduction = "harmony", dims = 1:25)
DimPlot(sce.all.harmony, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T)
DimPlot(sce.all.harmony, group.by = "tissue", reduction = "umap", label = T)
DimPlot(sce.all.harmony, group.by = "disease", reduction = "umap", label = T)
DimPlot(sce.all.harmony, group.by = "Individual", reduction = "umap", label = T)
#T-SNE
# sce.all <- RunTSNE(sce.all, dims = 1:25, check_duplicates = FALSE)
# DimPlot(sce.all, reduction = "tsne", label = T) 

save(sce.all.harmony, file = "GSE121893/GSE121893-sce.all.harmony.Rdata")  
# Annotation--------------------------------------------------------------------
sce.all.harmony <- RegroupIdents(sce.all.harmony, metadata = "RNA_snn_res.0.3")
ACM <- c('TNNT2', 'TTN', 'MYH6', 'MYH7') # 0
VCM <- c('MYH7', 'MYL2', 'FHL2') # 0
EC <- c('VWF', 'CDH5', 'THBD', 'TEK', 'ENG', 'PECAM1') # 1 5
FB <- c('FN1', 'VIM', 'COL1A1', 'COL1A2') # 2 ?3
M <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74') # 6
Myeloid	<- c('CD163', 'CD14', 'C1QA') # 6
P <- c('HCN1', 'HCN4', 'CACNA1D') # pacemaker cells (P cells) ?0
Peri <- c('RGS5', 'ABCC9', 'KCNJ8') # 3
SMC <- c('MYH11', 'ACTA2','CNN1', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') # 4
Meso <- c('MSLN', 'WT1', 'BNC1') # 8
Neural <- c('PLP1', 'NRXN1', 'NRXN3') # 
Adi <- c('GPAM', 'FASN', 'LEP') #
Lymph <- c('CD3E', 'IL7R', 'CD40LG') # 7
Mast <- c('KIT', 'CPA3') #
RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ') # 9

VlnPlot(sce.all.harmony, group.by = "RNA_snn_res.0.3", features = RBC, slot = "counts", log = TRUE)

sce.all.harmony@meta.data$cell_type <- 0
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 0] <- "CM"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 1] <- "EC"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 2] <- "FB"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 3] <- "Peri"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 4] <- "SMC"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 5] <- "EC"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 6] <- "Myeloid"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 7] <- "Lymphoid"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 8] <- "Meso"
sce.all.harmony@meta.data$cell_type[sce.all.harmony@meta.data$RNA_snn_res.0.3 == 9] <- "RBC"

p1 <-DimPlot(sce.all.harmony, reduction = "umap", 
        group.by = "RNA_snn_res.0.3",
        label = TRUE, pt.size = 0.5) + NoLegend()
p2 <-DimPlot(sce.all.harmony, reduction = "umap", 
             group.by = "cell_type",
             label = TRUE, pt.size = 0.5) + NoLegend()
p3 <-DimPlot(sce.all.harmony, reduction = "umap", 
             group.by = "Individual",
             label = TRUE, pt.size = 0.5) + NoLegend()
p1 + p2 + p3

save(sce.all.harmony, file = "GSE121893/GSE121893-sce.all.harmony.annot.Rdata") 

# Differential Analysis---------------------------------------------------------
FeaturePlot(sce.all.harmony, reduction = "umap", features = Targetgene)
VlnPlot(sce.all.harmony, features = Targetgene, slot = "counts", log = TRUE, ncol = 1)

CM <- subset(sce.all.harmony, subset = cell_type == "CM")

FAF2 <- CM@assays$RNA$counts[rownames(CM@assays$RNA$counts) == Targetgene]
group <- CM@meta.data$disease

mean(FAF2[group == "N"])
mean(FAF2[group != "N"])
wilcox.test(FAF2 ~ group)

dat2pseudo <- data.frame(
  ID = CM@meta.data$Individual,
  disease = CM@meta.data$disease,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>% 
  filter(!duplicated(ID)) %>% 
  mutate(Targetgene = NA)

table(CM@meta.data$orig.ident)
# C1  C2  D1  D2  D4  D5  N1 N10 N11 N12 N13 N14  N2  N3  N4  N5  N6  N7  N8  N9 
# 593   1 143 206 321 110 379 186 282 224 252 214 648 638 223 417 362 137  73 274 

CM <- subset(CM, subset = Individual != "C2")

for (i in unique(CM@meta.data$Individual)) {
  tmp <- subset(x = CM, subset = Individual == i) 
  dat2pseudo$Targetgene[dat2pseudo$ID == i] <- 
    tmp@assays$RNA$counts[rownames(tmp@assays$RNA$counts) == Targetgene] %>% sum
  print(i)
}

dat2pseudo

dat2pseudo %>% 
  filter(ID != 'C2') %>% 
  group_by(disease) %>% 
  mutate(mean = mean(Targetgene)) %>% 
  distinct(mean)
wilcox.test(Targetgene ~ disease, data = dat2pseudo) 
