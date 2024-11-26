rm(list = ls());gc()

require(dplyr)
require(stringr)
require(Seurat)
require(SeuratObject)

# Reference---------------------------------------------------------------------
# https://lishensuo.github.io/posts/bioinfo/028%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--%E5%9F%BA%E4%BA%8E%E6%96%87%E7%8C%AE%E7%9A%84%E7%BB%86%E8%83%9E%E7%B1%BB%E5%9E%8B%E6%B3%A8%E9%87%8Amarker/
# http://117.50.127.228/CellMarker/
# http://www.360doc.com/content/12/0121/07/76149697_996050517.shtml
# https://cloud.tencent.com/developer/article/2401058

# * GSE109816 and GSE121893 [Normal & Heart Failure] ---------------------------
# https://www.nature.com/articles/s41556-019-0446-7#Sec9
# Samples from human participants
# This study was approved by the Ethics Committee of Fuwai Hospital, Chinese Academy of Medical Sciences and Peking Union Medical University, as well as by the Ministry of Science and Technology. Written informed consent for tissue donation, which clearly stated the purpose of our study, was obtained from all of the patients.
# 
# LV and/or LA tissue samples were collected from 14 previously healthy organ donors (12 males and 2 females). Right ventricular tissue samples were obtained from two additional previously healthy organ donors (male). All of the donors were deceased due to acute trauma or anoxia, and all of the hearts had a LV ejection fraction of greater than 50%, a cold ischaemia time of less than 6 h and were deemed to be not suitable for transplantation for technical or non-cardiac reasons. Atrial tissue samples (n = 9) were obtained from the human LA, whereas all of the ventricular tissue samples (n = 7) were dissected from the LV anterior wall close to the cardiac apex. The range of donor ages was 21–52 yr, with a median age of 45.5 yr. Information regarding the cells isolated from these tissue samples is provided in in Supplementary Table 1.
# 
# HF samples were collected from patients (6 males) undergoing heart transplantation. Atrial and ventricular tissue samples were dissected from the same regions as in normal hearts. Six cases of HF were included in our study, with two progressing from cHF and four progressing from dHF.
# 
# For heart samples from patients with HF who were treated with LVAD, two male patients with HF received continuous-flow LVAD implantation either as therapy or as a bridge to heart transplantation. The samples that were obtained before therapy were collected during the implantation procedure, which involved introducing an opening in the LV. Post-treatment samples were acquired during the removal of the device. The LV ejection fraction was less than 30% in both of the patients. One patient showed cardiac function recovery 166 d after LVAD implanting (LVAD1), whereas the other patient underwent cardiac transplantation 155 d after LVAD implantation owing to a lack of any signs of functional improvement (LVAD2). Ventricular tissue samples were collected during surgical operation of LVAD treatment.

dat1 <- data.table::fread("GSE121893_human_heart_sc_umi.csv.gz")
dat2 <- data.table::fread("GSE109816_normal_heart_umi_matrix.csv.gz")
dat  <- dat1 %>% 
  left_join(dat2, by = "V1")

dim(dat)
# [1] 25742 14928

## Created Seurat Object--------------------------------------------------------
raw.data <- dat %>%
  tibble::column_to_rownames("V1") %>%
  as.data.frame

metadata <- data.frame(
  barcode = colnames(raw.data),
  orig.ident = colnames(raw.data) %>% str_remove("[SC]*_") %>% str_remove("_[0-9]*_[0-9]*$"), 
  stringsAsFactors = FALSE
)

annot1 <- data.table::fread("C:/D/Labmates/陈乾乾/GSE121893/GSE121893_all_heart_cell_cluster_info.txt.gz")
annot2 <- data.table::fread("C:/D/Labmates/陈乾乾/GSE121893/GSE109816_normal_heart_cell_info.txt.gz")

unique(annot1$sample)
unique(annot2$Individual)
intersect(unique(annot1$sample),
          unique(annot2$Individual))

metadata <- metadata %>%
  left_join(annot1 %>%
              dplyr::select(ID, sample) %>% 
              mutate(barcode = ID),
            by = "barcode")

table(metadata$sample)
# C1   C2   D1   D2   D4   D5   N1  N10  N11  N12  N13  N14   N2   N3   N4   N5   N6   N7 
# 1162  278  788 1171  601  221  654  470  511  381  466  246  624  607  259  736 1051  196 
# N8   N9 
# 287  668 

metadata <- metadata %>% dplyr::select(barcode, sample)
colnames(metadata)[2] <- "orig.ident"

all(metadata$barcode == colnames(raw.data))
# [1] TRUE
rownames(metadata) <- metadata$barcode

# key point
# https://www.biostars.org/p/9568513/
# To avoid ERROR: 'VST.default' is not implemented yet
raw.data <- as.matrix(raw.data)
raw.data <- Matrix(raw.data, sparse = TRUE)

sce.all <- CreateSeuratObject(counts = raw.data, data = raw.data)
sce.all <- AddMetaData(object = sce.all, metadata = metadata)
## QC --------------------------------------------------------------------------
# 质控的参数主要有两个： 
# 1.每个细胞测到的unique feature数目（unique feature代表一个细胞检测到的基因的数目，可以根据数据的质量进行调整） 
# 2.每个细胞检测到的线粒体基因的比例，理论上线粒体基因组与核基因组相比，只占很小一部分。所以线粒体基因表达比例过高的细胞会被过滤。
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all[["percent.rp"]] <- PercentageFeatureSet(sce.all, pattern = "^RP[SL]")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.rp"), ncol = 3)
# nCount代表每个细胞测到所有基因的表达量之和，
# percent.mt代表测到的线粒体基因的比例。

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

sce.all <- subset(sce.all, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.rp < 10)  
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.rp"), ncol = 3)

## Normalization ---------------------------------------------------------------
sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize", scale.factor = 10000)
#鉴定细胞间表达量高变的基因（feature selection）
#这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
#我们使用默认参数，即“vst”方法选取2000个高变基因。
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sce.all), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce.all)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

## Classification---------------------------------------------------------------
#### 1) Dimension Reduction-----------------------------------------------------
#Scaling the data
all.genes <- rownames(sce.all)
# Scales and centers features in the dataset. 
# If variables are provided in vars.to.regress, they are individually regressed against each feature, and the resulting residuals are then scaled and centered.
sce.all <- ScaleData(sce.all, features = all.genes)
#Perform linear dimensional reduction
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
#Examine and visualize PCA results a few different ways
print(sce.all[["pca"]], dims = 1:10, nfeatures = 5)
VizDimLoadings(sce.all, dims = 1:2, reduction = "pca")
DimPlot(sce.all, reduction = "pca")
DimHeatmap(sce.all, dims = 1:4, cells = 500, balanced = TRUE, )
#### 2) Determine Dims----------------------------------------------------------
#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
# JackStraw和Elbow都可以决定数据的“维度”。但是Elbow比较直观，我们选择Elbow结果进行解读。可以看到，主成分（PC）7到10之间，数据的标准差基本不在下降。所以我们需要在7到10之间进行选择，为了尊重官网的建议，我们选取10，即前10个主成分用于细胞的分类。
sce.all <- JackStraw(sce.all, dims = 50, num.replicate = 100) ## at most 50 PCs
sce.all <- ScoreJackStraw(sce.all, dims = 1:50)
JackStrawPlot(sce.all, dims = 1:50)
ElbowPlot(sce.all, ndims = 50)

#### 3) Cluster--------------------------------------------------------------
# 选择不同的resolution值可以获得不同的cluster数目，值越大cluster数目越多，默认值是0.8.

sce.all <- FindNeighbors(sce.all, dims = 1:25)
sce.all <- FindClusters(sce.all, resolution = 0.5)

#Look at cluster IDs of the first 5 cells
head(Idents(sce.all), 5)

save(sce.all, file = "GSE121893(PCA_1-25_0.5).Rdata")
#### 4) Visualization-----------------------------------------------------------
load("GSE121893(PCA).Rdata")
# TSNE和UMAP两种方法经常被用于可视化细胞类别。
#UMAP
sce.all <- RunUMAP(sce.all, dims = 1:25)
# head(sce.all@reductions$umap@cell.embeddings) 
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p1 <- DimPlot(sce.all, reduction = "umap", label = T)
#T-SNE
sce.all <- RunTSNE(sce.all, dims = 1:25, check_duplicates = FALSE)
# head(sce.all@reductions$tsne@cell.embeddings)
p2 <- DimPlot(sce.all, reduction = "tsne", label = T)
# p1 + p2

## Find markers-----------------------------------------------------------------
# 利用 FindMarkers 命令，可以找到找到各个细胞类型中与其他类别的差异表达基因，作为该细胞类型的生物学标记基因。
# 其中ident.1参数设置待分析的细胞类别，min.pct表示该基因表达数目占该类细胞总数的比例

#find all markers of cluster 1
# http://www.360doc.com/content/12/0121/07/76149697_996050517.shtml
# 参考曾老师的讲解，FindMarkers这里实际上是对第一类细胞和其余所有细胞做wilcox检验
# 返回差异基因。默认使用高变异基因，如果想使用全部细胞做差异检验：
# markers <- FindMarkers(sce.sub.Endo ,  ident.1="HE", ident.2="CI"
# assay = 'RNA',slot = 'counts',
# logfc.threshold =0,min.pct = 0 )
cluster3.markers <- FindMarkers(sce.all, ident.1 = "Unknown", min.pct = 0.25)
# head(cluster1.markers, n = 5)
#利用 DoHeatmap 命令可以可视化marker基因的表达
sce.all.markers <- FindAllMarkers(sce.all, only.pos = TRUE, min.pct = 0.25)
# ?FindMarkers
top3 <- sce.all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DoHeatmap(sce.all, features = top3$gene) + NoLegend()

## Exploration and Annotation---------------------------------------------------
VlnPlot(sce.all, features = c("MS4A1", "CD79A"))
VlnPlot(sce.all, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(sce.all, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

if(TRUE){
  ACM <- c('TNNT2', 'TTN', 'MYH6', 'MYH7') # 0 1 5
  VCM <- c('MYH7', 'MYL2', 'FHL2')
  EC <- c('VWF', 'CDH5', 'THBD', 'TEK', 'ENG', 'PECAM1') # 2 8
  FB <- c('FN1', 'VIM', 'COLA1', 'COL1A2') # 3
  M <- c('CD163', 'S100A8', 'CSF1R', 'C5AR1', 'CD74') # 9
  Myeloid	<- c('CD163', 'CD14', 'C1QA') # 9
  P <- c('HCN1', 'HCN4', 'CACNA1D') # pacemaker cells (P cells)
  Peri <- c('RGS5', 'ABCC9', 'KCNJ8')
  SMC <- c('MYH11', 'ACTA2','CNN1', 'TAGLN', 'MYOCD', 'CALD1', 'MYLK') # 6
  Meso <- c('MSLN', 'WT1', 'BNC1') # 14
  Neural <- c('PLP1', 'NRXN1', 'NRXN3')
  Adi <- c('GPAM', 'FASN', 'LEP')
  Lymph <- c('CD3E', 'IL7R', 'CD40LG') # 10
  Mast <- c('KIT', 'CPA3')
  RBC <- c('HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ')
  
  VlnPlot(sce.all, features = RBC, slot = "counts", log = TRUE)
}

new.cluster.ids <- c("CM", "CM", "EC", "FB", "Peri", "CM", 
                     "SMC", "Unknown", "EC", "M", "Lymph", 
                     "CM", "RBC", "CM", "Meso")
names(new.cluster.ids) <- levels(sce.all)
sce.all <- RenameIdents(sce.all, new.cluster.ids)
DimPlot(sce.all, reduction = "tsne", 
        # group.by = "orig.ident", 
        label = TRUE, pt.size = 0.5) + NoLegend()

unknown.markers <- FindMarkers(sce.all, ident.1 = "Unknown", min.pct = 0.25)

# Idents(sce.all) <- sce.all@meta.data$RNA_snn_res.0.5
DimPlot(sce.all, reduction = "tsne", 
        # group.by = "orig.ident", 
        label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(sce.all, reduction = "tsne", features = c("FAF2"))
VlnPlot(sce.all, features = c("FAF2"), slot = "counts", log = TRUE, ncol = 2)

sce.all.annot <- subset(x = sce.all, idents = c("CM", "EC", "FB", "Peri","SMC", "M", "Lymph", "Meso"))

DimPlot(sce.all, reduction = "tsne", 
        # group.by = "orig.ident", 
        label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(sce.all.annot, reduction = "tsne", features = c("CD163"))

VlnPlot(sce.all.annot, features = c("CD163"), slot = "counts", log = TRUE, ncol = 1)

# Note:
# 会影响分群和UMAP可视化的参数：
# 从这次的结果看，增加高变基因数目会把T细胞和B细胞分开。而PCA维数的增加没有使T细胞和B细胞分开，但是把T细胞分成了两个部分。
# 之前在这篇推文初探单细胞分析 — 标准化与降维聚类分群的理解提到过：高变基因个数越多，PCs的数目越多，保留的数据信息也就越多，更有可能引入噪音，运行速度也越慢。但是也不能太少，否则会丢失很多数据信息。
# 所以有可能是增加的高变基因使得T细胞和B细胞之间有更多的差异，可以分开。而PCA维数的增加把T细胞分成了两个部分，可能是由于有噪音引入。

# 不会影响分群，只会影响UMAP可视化的参数：
# UMAP参数中的n_neighbors，min_dist和dims这三个参数，是不会影响细胞的本质属性的，只是影响UMAP可视化的图。可以看到减小的min_dist和增加的dims会把T细胞和B细胞分开。n_neighbors的影响不大。
