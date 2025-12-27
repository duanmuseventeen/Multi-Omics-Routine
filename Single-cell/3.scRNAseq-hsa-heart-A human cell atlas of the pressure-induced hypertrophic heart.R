# A human cell atlas of the pressure-induced hypertrophic heart ---------------
# https://www.nature.com/articles/s44161-023-00260-8
# https://github.com/djhn75/DiseaseHeartCellAtlas

# Reference
# https://zhuanlan.zhihu.com/p/614993353

# Note:
# Even someone thinks we cannot use the Rdata directly, rather than we should 
# run cluster again by ourselves. Here, I do not reanalyze the data, because this
# analysis aims to validate our results from wet experiment.

# Load data
# Download from https://explore.data.humancellatlas.org/projects/902dc043-7091-445c-9442-d72e163b9879
sce.all <- readRDS("A human cell atlas of the pressure-induced hypertrophic heart/seurat_object_hca_as_harmonized_AS_SP_nuc_refined_cells.rds")

sce.all <- RegroupIdents(sce.all, metadata = "cell_type")
# Filter low quality samples
sce.all.annot <- subset(sce.all, subset = cell_type %in% c("CM", "PC", "SMC", "EC","Lympoid", "Myeloid", "Neuro", "FB"))

# Visualization
DimPlot(sce.all.annot, reduction = "umap",
        group.by = "cell_type",
        label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(sce.all.annot, features = c(Targetgene))
VlnPlot(object = sce.all.annot, features = Targetgene, slot = "counts", log = TRUE, group.by = "cell_type")

# Subset Cardial Myocytes
CM <- subset(sce.all, subset = cell_type %in% c("CM"))

# Differential Analysis for target gene 1
Targetgene <- CM@assays$RNA$counts[rownames(CM@assays$RNA$counts) == Targetgene]
group <- CM@meta.data$dataset

VlnPlot(object = CM, features = Targetgene, slot = "counts", log = TRUE, group.by = "dataset")
mean(FAF2[group == "hca"])
mean(FAF2[group != "hca"])
wilcox.test(FAF2 ~ group)
# Differential Analysis for target gene 2
Normal <- subset(CM, subset = dataset %in% c("hca"))
HF     <- subset(CM, subset = dataset %in% c("aort_sten"))

Normal.Target <- Normal@assays$RNA$counts[rownames(Normal@assays$RNA$counts) == Targetgene]
HF.Target     <- HF@assays$RNA$counts[rownames(HF@assays$RNA$counts) == Targetgene]

mean(Normal.Target)
mean(HF.Target)
wilcox.test(Normal.Target, HF.Target)

# Differential Analysis for all genes
CM <- RegroupIdents(CM, metadata = "dataset")
FindMarkers(CM, group.by = "dataset", ident.1 = 'hca', ident.2 = 'aort_sten', logfc.threshold = 0, min.pct = 0)

# Pseudobulk Differential Analysis
dat2pseudo <- data.frame(
  ID = CM@meta.data$orig.ident,
  disease = CM@meta.data$dataset,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>% 
  filter(!duplicated(ID)) %>% 
  mutate(Targetgene = NA)

for (i in unique(CM@meta.data$orig.ident)) {
  tmp <- subset(x = sce.all, subset = orig.ident == i) 
  dat2pseudo$Targetgene[dat2pseudo$ID == i] <- 
    tmp@assays$RNA$counts[rownames(tmp@assays$RNA$counts) == Targetgene] %>% sum
  print(i)
}

dat2pseudo %>% group_by(disease) %>% mutate(mean = mean(Targetgene)) %>% distinct(mean)
wilcox.test(Targetgene ~ disease, data = dat2pseudo)
