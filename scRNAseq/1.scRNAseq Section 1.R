# A Proximal-to-Distal Survey of Healthy Adult Human Small Intestine and Colon Epithelium by Single-Cell Transcriptomics
# PMID: 35176508
# https://cellxgene.cziscience.com/collections/64b24fda-6591-4ce1-89e7-33eb6c43ad7b
# https://cellxgene.cziscience.com/e/019c7af2-c827-4454-9970-44d5e39ce068.cxg/

require(Seurat)
require(SeuratObject)
require(dplyr)

gse185224 <- readRDS("GSE185224.rds")

Idents(gse185224) <- gse185224@meta.data$cell_type

gse185224 %>%
  group_by(cell_type) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Find markers
markers <- FindAllMarkers(gse185224, only.pos = TRUE)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# Umap or T-SNE
DimPlot(gse185224, dims = c(1, 2), reduction = "umap", label = TRUE, pt.size = 0.5)

# Violin plot
VlnPlot(gse185224, features = c("ENSG00000111640","ENSG00000075624"))

# Feature plot
FeaturePlot(object = gse185224, features = 'PC_1')
FeaturePlot(object = gse185224, c("ENSG00000111640","ENSG00000075624")) # GAPDH ACTB
